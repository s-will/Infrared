#ifndef INFRARED_CLUSTER_TREE_HPP
#define INFRARED_CLUSTER_TREE_HPP

/*
 * InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
 * (C) Sebastian Will, 2018
 *
 * This file is part of the InfraRed source code.
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 */

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "constraint_network.hpp"

namespace ired {


    /**
     * @brief A tree of clusters (=variables, functions, constraints)
     *
     * Supports evaluation and sampling.
     *
     */
    template<class FunValue=double>
    class ClusterTree {

    public:
        using constraint_network_t = ConstraintNetwork<FunValue>;

        //! evaluation policy type
        using evaluation_policy_t = typename constraint_network_t::evaluation_policy_t;


        using var_idx_t = typename constraint_network_t::var_idx_t;
        using cluster_t = typename constraint_network_t::cluster_t;
        using assignment_t = typename constraint_network_t::assignment_t;
        using fun_value_t = typename constraint_network_t::fun_value_t;
        using constraint_t = typename constraint_network_t::constraint_t;
        using function_t = typename constraint_network_t::function_t;

        using message_t = MaterializedFunction<fun_value_t,mapS>;

        struct vertex_info_t {
            cluster_t cluster;
        };

        struct edge_info_t {
            Function<fun_value_t> *message;
        };

        using tree_t = boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, vertex_info_t, edge_info_t >;

        using vertex_descriptor = typename boost::graph_traits<tree_t>::vertex_descriptor;

        /**
         * @brief construct empty
         */
        ClusterTree() {
        }

        /**
         * @brief Construct with domains
         */
        ClusterTree(std::vector<int> &domsizes)
            : cn_(domsizes) {
        };

        /**
         * @brief Construct with uniform domains
         */
        ClusterTree(int num_vars, int domsize)
            : cn_(num_vars, domsize) {
        };


        //! @brief read access to constraint network
        const auto & constraint_network() const {
            return cn_;
        }

        //
        // enum class TDMethod {
        //     tdlib
        // };
        //
        //  /** @brief construct from constraint network by some tree decomposition method
        //  */
        // ClusterTree(constraint_network_t &n, TDMethod td_method): cn_(cn) {
        // }


        // for 'manual' construction
        auto
        add_root_cluster(const std::vector<var_idx_t> &vars) {
            auto node = boost::add_vertex(tree_);
            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        // auto
        // make_root_cluster_python(const boost::python::list &vars) {
        //     return make_root_cluster(python_extract_list<var_idx_t>(vars));
        // }

        auto
        add_child_cluster( vertex_descriptor parent, const std::vector<int> &vars) {
            auto node = boost::add_vertex(tree_);
            boost::add_edge(parent, node, tree_);

            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        // auto
        // make_child_cluster_python( vertex_descriptor parent, const boost::python::list &vars) {
        //     return make_child_cluster(parent, python_extract_list<var_idx_t>(vars));
        // }

        //! @brief add variable to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_variable( vertex_descriptor node, var_idx_t var ) {
            tree_[node].cluster.add_variable( var );
        }

        //! @brief add constraint to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_constraint( vertex_descriptor node, const std::shared_ptr<constraint_t> &x ) {
            tree_[node].cluster.add_constraint( cn_.add_constraint(x) );
        }

        //! @brief add function to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_function( vertex_descriptor node, const std::shared_ptr<function_t> &x ) {
            tree_[node].cluster.add_function( cn_.add_function(x) );
        }

        // run the dynamic programming evaluation
        void
        evaluate();

        // generate one sample, using the DP tables
        //
        // this is a real 'sample' only if partition functions are
        // computed due to the evaluation policy!
        auto
        sample();

    private:
        constraint_network_t cn_;
        tree_t tree_;

        bool evaluated_ = false;
        bool single_empty_rooted_ = false;

        vertex_descriptor root_;

        // insert pseudo root to connect trees of the forest (unless already done)
        // @returns new root
        auto
        single_empty_root();

        struct evaluate_finish_edge {
            using event_filter = boost::on_finish_edge;

            evaluate_finish_edge(constraint_network_t &cn, tree_t &tree)
                : cn_(cn), tree_(tree) {
            }

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {

                const auto &parent = graph[ source(e, graph) ].cluster;
                const auto &child =  graph[ target(e, graph) ].cluster;

                // compute the message from cluster child to cluster parent

                auto sep  = child.sep_vars(parent);
                auto diff = child.diff_vars(parent);

                // concat sep + diff
                auto sep_diff = sep;
                sep_diff.insert(sep_diff.end(),diff.begin(),diff.end());

                auto message = std::make_unique<message_t>(sep, cn_);

                auto a = Assignment(cn_.num_vars());
                
                auto it = a.make_iterator
                    (sep_diff,
                     cn_,
                     child.constraints(),
                     child.functions(),
                     //evaluate 0-ary functions
                     a.eval_determined(child.functions(), evaluation_policy_t()) 
                     );

                fun_value_t x = evaluation_policy_t::zero();

                it.register_finish_stage2_hook
                    (sep.size(),
                     [&message,&a,&x] () { 
                        message->set(a, x);
                        x = evaluation_policy_t::zero();
                    } );

                for(; ! it.finished() ; ++it ) {
                    x = evaluation_policy_t::plus( x, it.value() );
                }

                // register message in cn, such that it persists!
                auto msg = cn_.add_function(std::move(message));
                // then, register in cluster parent
                tree_[ source(e, graph) ].cluster.add_function(msg);
                
                // ... and as edge property
                tree_[e].message = msg;
            }

        private:
            constraint_network_t &cn_;
            tree_t &tree_;
        };

        struct sample_examine_edge {
            using event_filter = boost::on_examine_edge;

            sample_examine_edge(constraint_network_t &cn,assignment_t &a) 
                : cn_(cn), a_(a) {}

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {
                const auto &parent = graph[ source(e, graph) ].cluster;
                const auto &child =  graph[ target(e, graph) ].cluster;

                auto diff = child.diff_vars(parent);

                const auto &message = graph[e].message;
                
                double r = rand()/(RAND_MAX+1.0) * (*message)(a_);

                assert ( a_.eval_determined(child.constraints(), StdEvaluationPolicy<bool>()) );
                
                a_.set_undet(diff);

                auto x = evaluation_policy_t::zero();
                auto it = a_.make_iterator(diff, cn_,
                                           child.constraints(),
                                           child.functions(),
                                           a_.eval_determined(child.functions(), evaluation_policy_t())
                                           );
                for( ; ! it.finished(); ++it ) {
                    x = evaluation_policy_t::plus( x, it.value() );
                    
                    if ( x > r ) {
                        break;
                    }
                }
            }

        private:
            constraint_network_t &cn_;
            assignment_t &a_;
        };

    };

    template<class ConstraintNetwork>
    auto
    ClusterTree<ConstraintNetwork>::single_empty_root() {
        if (single_empty_rooted_) {return root_;}
        
        // find all root nodes of the tree
        auto index = boost::get(boost::vertex_index, tree_);

        std::set<int> no_roots;
        std::vector<vertex_descriptor> old_roots;

        for (auto it = edges(tree_).first; it != edges(tree_).second; ++it) {
            no_roots.insert(index[target(*it,tree_)]);
        }

        for (auto it = vertices(tree_).first; it != vertices(tree_).second; ++it) {
            if ( no_roots.find(index[*it]) == no_roots.end() ) {
                old_roots.push_back(*it);
            }
        }
        
        if ( old_roots.size()==1 && tree_[ old_roots[0] ].cluster.empty() ) {
            root_ = old_roots[0];
        } else {
            // insert new root and point to all old roots
            root_ = boost::add_vertex(tree_);
            
            for (auto old_root : old_roots) {
                boost::add_edge(root_, old_root , tree_);
            }
        }

        single_empty_rooted_ = true;
        
        return root_;
    }

    template<class ConstraintNetwork>
    void
    ClusterTree<ConstraintNetwork>
    ::evaluate() {

        auto root = single_empty_root();

        auto cte_visitor = boost::make_dfs_visitor(evaluate_finish_edge(cn_,tree_));

        boost::depth_first_search(tree_, visitor(cte_visitor).root_vertex(root));

        evaluated_ = true;
    }

    template<class ConstraintNetwork>
    auto
    ClusterTree<ConstraintNetwork>::sample() {

        if (!evaluated_) {
            evaluate();
        }

        auto a = assignment_t(cn_.num_vars());

        auto sample_visitor = boost::make_dfs_visitor(sample_examine_edge(cn_,a));

        boost::depth_first_search(tree_, visitor(sample_visitor).root_vertex(single_empty_root()));

        return a;
    }
}

#endif
