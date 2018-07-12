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
    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy<FunValue>>
    class ClusterTree {

    public:
        using constraint_network_t = ConstraintNetwork<FunValue>;

        using cluster_t = typename constraint_network_t::cluster_t;
        using assignment_t = typename constraint_network_t::assignment_t;
        using fun_value_t = typename constraint_network_t::fun_value_t;
        using constraint_t = typename constraint_network_t::constraint_t;
        using function_t = typename constraint_network_t::function_t;


        using message_t = MaterializedFunction<fun_value_t>;


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
            : cn_(num_vars,domsize) {
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
        make_root_cluster(const std::vector<int> &vars) {
            auto node = boost::add_vertex(tree_);
            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        auto
        make_child_cluster( vertex_descriptor parent, const std::vector<int> &vars) {
            auto node = boost::add_vertex(tree_);
            boost::add_edge(tree_, parent, node);

            tree_[node].cluster = cluster_t(vars);
            return node;
        }

        //! @brief add constraint to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_constraint( vertex_descriptor node, const constraint_t &c ) {
            tree_[node].cluster.add_constraint( cn_.add_constraint(c) );
        }

        //! @brief add function to cluster
        //!
        //! synchronized between cluster and cn
        auto
        add_function( vertex_descriptor node, const function_t &c ) {
            tree_[node].cluster.add_function( cn_.add_function(c) );
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
        bool single_rooted_ = false;

        vertex_descriptor root_;

        // insert pseudo root to connect trees of the forest (if necessary)
        // @returns (potentially new) root
        auto
        single_root();

        struct evaluate_finish_edge {
            using event_filter = boost::on_finish_edge;

            evaluate_finish_edge(constraint_network_t &cn, tree_t &tree)
                : cn_(cn), tree_(tree) {
            }

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {

                auto parent = graph[ source(e, graph) ].cluster;
                auto child =  graph[ target(e, graph) ].cluster;


                // compute the message from cluster child to cluster parent

                auto sep  = child.sep_vars(parent);
                auto diff = child.diff_vars(parent);

                auto message = message_t(sep,cn_.domsizes());

                auto a = Assignment(cn_.num_vars());

                auto it = a.make_iterator(sep, cn_.domsizes(), child.constraints());

                for( ++it ; ! it.finished() ; ++it ) {
                    fun_value_t x = EvaluationPolicy::zero();

                    auto it2 = a.make_iterator(diff, cn_.domsizes(), child.constraints());

                    for( ++it2; ! it2.finished(); ++it2) {

                        fun_value_t p = EvaluationPolicy::one();
                        for ( auto f : child.functions() ) {
                            p = EvaluationPolicy::multiplies(p, (*f)(a) );
                        }
                        x = EvaluationPolicy::plus(x, p);
                    }

                    message.set(a, x);

                    a.reset(diff);
                }

                // register message in cn, such that it persists!
                auto msg = cn_.add_function(message);
                // then, register in cluster parent
                parent.add_function(msg);
                // ... and as edge property
                tree_[e].message = msg;
            }

        private:
            constraint_network_t &cn_;
            tree_t &tree_;
        };

        struct sample_examine_edge {
            using event_filter = boost::on_examine_edge;

            sample_examine_edge(constraint_network_t &cn,assignment_t &a) : cn_(cn), a_(a) {}

            template <class Edge, class Graph>
            void
            operator ()(Edge e, Graph &graph) {
                auto parent = graph[ source(e, graph) ].cluster;
                auto child =  graph[ target(e, graph) ].cluster;

                auto diff = child.diff_vars(parent);

                auto message = graph[e].message;

                double r = rand()/(RAND_MAX+1.0) * (*message)(a_);
                
                auto x = EvaluationPolicy::zero();
                for( auto it = ++ a_.make_iterator(diff, cn_.domsizes(), child.constraints()); ! it.finished(); ++it ) {
                    fun_value_t p = EvaluationPolicy::one();
                    for ( auto f : child.functions() ) {
                        p = EvaluationPolicy::multiplies(p, (*f)(a_) );
                    }

                    x = EvaluationPolicy::plus(x, p);

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

    template<class ConstraintNetwork, class EvaluationPolicy>
    auto
    ClusterTree<ConstraintNetwork, EvaluationPolicy>::single_root() {
        if (single_rooted_) {return root_;}

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

        // insert new root and point to all old roots
        auto new_root = boost::add_vertex(tree_);

        for (auto old_root : old_roots) {
            boost::add_edge(new_root, old_root , tree_);
        }

        single_rooted_ = true;
        root_ = new_root;
        return root_;
    }

    template<class ConstraintNetwork, class EvaluationPolicy>
    void
    ClusterTree<ConstraintNetwork,EvaluationPolicy>
    ::evaluate() {

        auto root = single_root();

        auto cte_visitor = boost::make_dfs_visitor(evaluate_finish_edge(cn_,tree_));

        boost::depth_first_search(tree_, visitor(cte_visitor).root_vertex(root));

        evaluated_ = true;
    }


    template<class ConstraintNetwork, class EvaluationPolicy>
    auto
    ClusterTree<ConstraintNetwork,EvaluationPolicy>::sample() {

        if (!evaluated_) {
            evaluate();
        }

        auto a = assignment_t(cn_.num_vars());

        auto sample_visitor = boost::make_dfs_visitor(sample_examine_edge(cn_,a));

        boost::depth_first_search(tree_, visitor(sample_visitor).root_vertex(single_root()));

        return a;
    }
}

#endif
