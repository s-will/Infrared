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

/**
 * @file
 *
 * @brief Defines the cluster tree.
 */

#include <set>
#include <optional>

#include "graph.hpp"
#include "feature_network.hpp"
#include "subtree_weights.hpp"

namespace ired {
    /**
     * @brief A tree of clusters (=variables, functions, constraints)
     *
     * A cluster tree belongs to a constraint network and holds
     * it. (The latter simplifies the usage in Infrared; however this
     * means, we do not support several cluster trees for the same
     * constraint network.)
     *
     * The cluster tree supports evaluation and sampling. For this
     * purpose, the cluster tree must satisfy the cluster tree
     * properties
     *
     * 1) For each variable in the CN, there is one bag that contains
     * the variable.
     *
     * 2) For each variable, the bags that contain this variable form
     * a subtree.
     *
     * 3) For each function, there is exactly one bag that contains
     * the function and its variables.
     *
     * 4) For each constraint, there is at least one bag that contains
     * the function and its variables. Constraints are only assigned to bags
     * that contain all of their variables.
     *
     * Ensuring those properties is the job of the user of this class!
     */
    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy<FunValue>>
    class ClusterTree {

    public:
        using feature_network_t = FeatureNetwork<FunValue, EvaluationPolicy>;

        //! evaluation policy type
        using evaluation_policy_t = typename feature_network_t::evaluation_policy_t;

        using var_idx_t = typename feature_network_t::var_idx_t;
        using cluster_t = typename feature_network_t::cluster_t;
        using assignment_t = typename feature_network_t::assignment_t;
        using fun_value_t = typename feature_network_t::fun_value_t;
        using constraint_t = typename feature_network_t::constraint_t;
        using function_t = typename feature_network_t::function_t;

        using message_t = MaterializedFunction<fun_value_t>;

        //! @brief information at a vertex (=cluster/bag) of the tree
        struct vertex_info_t {
            vertex_info_t()
                : cluster() {}
            vertex_info_t(const cluster_t &cluster)
                : cluster(cluster) {}
            cluster_t cluster;
        };

        //! @brief information at an edge of the tree
        struct edge_info_t {
            edge_info_t(function_t *message=nullptr)
                : message(message) {}
            function_t *message;
        };

        using tree_t = graph::adjacency_list<vertex_info_t, edge_info_t>;

        //! @brief type of identifiers of vertices (typically 'long int')
        using vertex_descriptor_t = typename tree_t::vertex_descriptor_t;

        using weighted_function_t = WeightedFunction<fun_value_t>;
        using weights_tree_t = WeightsTree<fun_value_t, evaluation_policy_t>;
        using cluster_info_t = typename weights_tree_t::cluster_info_t;

        /**
         * @brief construct empty
         */
        ClusterTree() {
        }

        /**
         * @brief Construct with variable domains
         *
         * @param domsizes vector of domain sizes for each variable;
         * its length specifies the number of variables
         */
        explicit
        ClusterTree(const FiniteDomainVector &domains)
            : fn_( domains ) {
        }

        /**
         * @brief construct from vector of upper bounds
         */
        explicit
        ClusterTree(std::vector<int> domsizes)
            : fn_( domains_from_domsizes( domsizes ) ) {
        }

        /**
         * @brief Construct with uniform domains
         *
         * @param num_vars the number of variables in the underlying constraint network
         * @param domsize uniform domain size of each variable
         */
        ClusterTree(int num_vars, const FiniteDomain &domain)
            : fn_(num_vars, domain) {
        }

        /**
         * @brief Construct with uniform domains
         *
         * @param num_vars the number of variables in the underlying constraint network
         * @param domsize uniform domain size of each variable
         */
        ClusterTree(int num_vars, int domsize)
            : fn_( num_vars, FiniteDomain(0, domsize-1) ) {
        }

        ~ClusterTree() {}

        //! @brief read access to constraint network
        const auto & feature_network() const {
            return fn_;
        }

        /**
         * @brief add new root cluster to the tree
         *
         * @param vars variables of the new cluster
         * @return bag/cluster id
         *
         * Typically used when constructing the cluster tree. Typically, one
         * adds functions and constraints to the new cluster that is
         * created by this method (using the returned bag id).
         */
        auto
        add_root_cluster(const std::vector<var_idx_t> &vars) {
            return tree_.add_vertex(cluster_t(vars));
        }

        /**
         * @brief add new child cluster to the tree
         *
         * @param parent the parent of this child
         * @param vars variables of the new cluster
         * @return bag/cluster id
         *
         * Typically used when constructing the cluster tree. Typically, one
         * adds functions and constraints to the new cluster that is
         * created by this method (using the returned bag id).
         */
        auto
        add_child_cluster( vertex_descriptor_t parent, const std::vector<int> &vars) {
            auto child = tree_.add_vertex(cluster_t(vars));
            tree_.add_edge(parent, child, edge_info_t());
            return child;
        }

        /**
         * @brief add new constraint to cluster
         *
         * @param id identifier of cluster
         * @param x (shared pointer to) constraint
         *
         * Typically used when constructing the cluster tree.
         */
        void
        add_constraint( vertex_descriptor_t node, const std::shared_ptr<constraint_t> &x ) {
            tree_[node].cluster.add_constraint( fn_.add_constraint(x) );
        }

        /**
         * @brief add new function to cluster
         *
         * @param id identifier of cluster
         * @param x (shared pointer to) function
         *
         * Typically used when constructing the cluster tree.
         */
        void
        add_function( vertex_descriptor_t node, const std::shared_ptr<function_t> &x ) {
            tree_[node].cluster.add_function( fn_.add_function(x) );
        }

        /**
         * @brief Evaluate the cluster tree (by DP)
         *
         * @return evaluation result
         *
         * Call this once before generating tracebacks with traceback()
         *
         * @note multiple calls will not rerun the evaluation algorithm,
         * but return the stored evaluation result
         */
        auto
        evaluate();

        /**
         * @brief Check consistency
         *
         * @return whether the ct is consistent
         *
         * @note Implicitely evaluates the cluster tree by calling @see
         * evaluate().  Consequently, after this method returned true,
         * valid tracebacks can be produced by @see traceback().
         */
        bool is_consistent() {
            return evaluate() != evaluation_policy_t::zero();
        }

        /**
         * @brief Generate a traceback
         *
         * @return assignment obtained by traceback
         *
         * @note The nature of this traceback depends on the evaluation
         * policy!
         *
         * Generates one traceback assignment. Requires that the cluster tree
         * was evaluated and is consistent.
         * Arbitrarily many tracebacks can be generated after a single
         * evaluation (which is of course mostly reasonable for non-det tb,
         * like stochastic traceback)
         *
         * Results are NOT defined if the cluster tree was never evaluated (due to a call
         * to evaluate() or is_consistent()) or is not consistent. Calling
         * traceback on an inconsistent tree may even result in termination.
         *
         * On evaluated and consistent trees, stochastic traceback, samples from a distribution
         * controlled by the functions and constraints. When the functions
         * return Boltzmann weights, assignments are traced back from a
         * Boltzmann distribution.
         */
        auto traceback();

        /**
         * @brief Generate a restricted traceback
         *
         * @param variables list of variables to which trace is restricted
         * @param assignment assignment of the other variables
         *
         * @return assignment obtained by restricted traceback
         *
         * @see traceback
         *
         */
        auto restricted_traceback(const std::set<var_idx_t> &variables,
            const assignment_t &assignment);

        /**
         * @brief Generate a (possibly non-redundant) traceback
         *
         * @param non_redundant specifies whether the samples should be unique
         * @param non_redundant_mode type of non-redundant algorithm
         *
         * @return assignment obtained by (possibly non-redundant) traceback
         *
         * If non_redundant is false, operates in the same way as @see traceback.
         * 
         * If non_redundant_mode == "rejection", stores all generated assignments. 
         * If an assignment is already present in the list throws an exception.
         * 
         * If non_redundant_mode == "default", uses subtree weights 
         * to ensure the generated assignment are unique. In case 
         * there are no more possible assignments throws an exception.
         */
        auto traceback_new(bool non_redundant, std::string non_redundant_mode);

    private:
        feature_network_t fn_;
        tree_t tree_;

        bool evaluated_ = false;
        fun_value_t evaluation_result_;
        bool single_empty_rooted_ = false;

        vertex_descriptor_t root_;

        std::optional<std::set<var_idx_t>> trace_variables;

        std::vector<assignment_t> forbidden;
        weights_tree_t weights_tree_;
        
        // insert pseudo root to connect trees of the forest (unless already done)
        // @returns new root
        auto
        single_empty_root();

        auto domains_from_domsizes( std::vector<int> &domsizes ) {
            auto domains = FiniteDomainVector();
            for (auto x: domsizes) {
                domains.push_back( FiniteDomain( 0, x - 1 ) );
            }
            return domains;
        }

        void
        dfs_evaluate(vertex_descriptor_t v) {
            // notes: this terminates if there are no children
            // as well, we can know, that tree_ is a tree, without any
            // cycles

            for(auto &e: tree_.adj_edges(v)) {
                // recursively evaluate subtree at e.target()
                dfs_evaluate(e.target());

                // then, compute the message for the current edge
                // from cluster child to cluster parent

                const auto &parent = tree_[ e.source() ].cluster;
                const auto &child =  tree_[ e.target() ].cluster;

                auto sep  = child.sep_vars(parent);
                auto diff = child.diff_vars(parent);

                // concat sep + diff
                auto sep_diff = sep;
                sep_diff.insert(sep_diff.end(),diff.begin(),diff.end());

                auto message = std::make_unique<message_t>(sep, fn_);

                auto a = assignment_t(fn_.domains());

                auto it = a.make_iterator
                    (sep_diff,
                     fn_,
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
                auto msg = fn_.add_function(std::move(message));
                // then, register in cluster parent
                tree_[ e.source() ].cluster.add_function(msg);

                // ... and as edge property
                e.message = msg;
            }
        }

        /**
         * @brief trace back by dfs
         * @param v trace back from this node in tree_
         * @param a [inout] the assignment that determines the variables above
         * v; as well used to return the variables from trace back in the
         * subtree of v
         */
        void
        dfs_traceback(vertex_descriptor_t v, assignment_t &a) {

            for(const auto &e: tree_.adj_edges(v)) {

                const auto &parent = tree_[ e.source() ].cluster;
                const auto &child =  tree_[ e.target() ].cluster;

                auto diff = child.diff_vars(parent);
                bool diff_changed = false;

                if (trace_variables) { // if the trace is restricted,
                                       // then remove other vars from diff
                    auto new_end = remove_if(diff.begin(), diff.end(),
                            [&](var_idx_t i) {
                                return trace_variables->find(i) == trace_variables->end();
                            });
                    diff_changed = new_end != diff.end();
                    diff.erase(new_end, diff.end());
                }

                if (! diff.empty()) {

                    a.set_undet(diff);

                    assert(a.eval_determined(child.constraints(),
                                             StdEvaluationPolicy<bool>()));

                    auto it = a.make_iterator(diff, fn_,
                                              child.constraints(),
                                              child.functions(),
                                              a.eval_determined(child.functions(),
                                                                evaluation_policy_t())
                                              );

                    // -----
                    // determine value for selector
                    auto value = evaluation_policy_t::zero();
                    if (diff_changed) {
                        // recompute value for selector
                        for( ; ! it.finished(); ++it ) {
                            value = evaluation_policy_t::plus( value, it.value() );
                        }
                        it.reset();
                    } else {
                        // get value by evaluating message:
                        // evaluate message at partial assignment a;
                        value = (*e.message)(a);
                    }
                    // initialize the Selector
                    auto selector = typename evaluation_policy_t::selector(value);

                    // enumerate all combinations of values assigned to diff variables;
                    // recompute SUM_a(PROD_f(f(a))), where a runs over these
                    // combinations
                    auto x = evaluation_policy_t::zero();
                    for( ; ! it.finished(); ++it ) {
                        x = evaluation_policy_t::plus( x, it.value() );

                        if (selector.select(x)) {
                            break;
                        }
                    }
                }
                // trace back from target
                dfs_traceback(e.target(), a);
            }
        }

        /**
         * @brief trace back by dfs
         * @param v trace back from this node in tree_
         * @param a [inout] the assignment that determines the variables above
         * v; as well used to return the variables from trace back in the
         * subtree of v
         * @param non_redundant specifies whether subtree weights should be added 
         * to ensure non-redundancy
         */
        void
        dfs_traceback_new(vertex_descriptor_t v, assignment_t &a, bool non_redundant) {

            for(const auto &e: tree_.adj_edges(v)) {

                const auto &parent = tree_[ e.source() ].cluster;
                const auto &child =  tree_[ e.target() ].cluster;

                auto diff = child.diff_vars(parent);
                bool diff_changed = false;

                if (trace_variables) { // if the trace is restricted,
                                       // then remove other vars from diff
                    auto new_end = remove_if(diff.begin(), diff.end(),
                            [&](var_idx_t i) {
                                return trace_variables->find(i) == trace_variables->end();
                            });
                    diff_changed = new_end != diff.end();
                    diff.erase(new_end, diff.end());
                }

                if (! diff.empty()) {

                    a.set_undet(diff);

                    assert(a.eval_determined(child.constraints(),
                                             StdEvaluationPolicy<bool>()));

                    // add weights to the messages from the children nodes
                    if (non_redundant) {
                        size_t num_messages = tree_.adj_edges(e.target()).size();
                        size_t total_num_funs = child.functions().size();
                        for (size_t i=total_num_funs-num_messages; i<total_num_funs; ++i) {
                            fun_value_t Z_c = weights_tree_.get_weight(e.target(), a);
                            function_t *weighted_function = new weighted_function_t(child.functions()[i]->vars(), child.functions()[i], Z_c);
                            child.modify_function(i, weighted_function);
                        }
                    }

                    auto it = a.make_iterator(diff, fn_,
                                              child.constraints(),
                                              child.functions(),
                                              a.eval_determined(child.functions(),
                                                                evaluation_policy_t())
                                              );
                    
                    // -----
                    // determine value for selector
                    auto value = evaluation_policy_t::zero();
                    if (diff_changed) {
                        // recompute value for selector
                        for( ; ! it.finished(); ++it ) {
                            value = evaluation_policy_t::plus( value, it.value() );
                        }
                        it.reset();
                    } else {
                        // get value by evaluating message:
                        // evaluate message at partial assignment a;
                        value = (*e.message)(a);
                    }
                    // add weight to the value for selector
                    if (non_redundant) {
                        fun_value_t Z = weights_tree_.get_weight(e.source(), a);
                        value = evaluation_policy_t::plus(value, -Z);
                        if (value==evaluation_policy_t::zero()) {
                            // if the value becomes zero, there are no more new samples
                            throw std::exception();
                        }
                    }
                    

                    // initialize the Selector
                    auto selector = typename evaluation_policy_t::selector(value);

                    // enumerate all combinations of values assigned to diff variables;
                    // recompute SUM_a(PROD_f(f(a))), where a runs over these
                    // combinations
                    auto x = evaluation_policy_t::zero();
                    for( ; ! it.finished(); ++it ) {
                        x = evaluation_policy_t::plus( x, it.value() );

                        if (selector.select(x)) {
                            break;
                        }
                    }
                }
                // trace back from target
                dfs_traceback_new(e.target(), a, non_redundant);
            }
        }

        /**
         * @brief Traverse the cluster tree and store its nodes in order
         * 
         * @param n current node in the cluster tree
         * @param nodes list of collected node information 
         * (cluster and children) for the previous nodes
         * @param p possition of the current node in df order
         * @param preorder list of possitions corresponding 
         * to the vector describtors from the cluster tree
         */
        void
        df_traversal(vertex_descriptor_t n, std::vector<cluster_info_t> &nodes,
                     vertex_descriptor_t &p, std::vector<vertex_descriptor_t> &preorder) {
            std::vector<vertex_descriptor_t> children;
            for(const auto &e: tree_.adj_edges(n)) {
                children.push_back(e.target());
            }

            cluster_info_t item(tree_[n].cluster, children);
            nodes.push_back(item);

            preorder[n] = p;
            p++;
            
            for(const auto &e: tree_.adj_edges(n)) {
                df_traversal(e.target(), nodes, p, preorder);
            }
        }

    }; // end class ClusterTree

    // return single empty cluster that roots the tree; if such a
    // cluster exists or was generated before, simply return it;
    // otherwise, construct a new empty cluster, connect it as parent
    // to all existing roots, and return it.
    template<class FunValue, class EvaluationPolicy>
    auto
    ClusterTree<FunValue,EvaluationPolicy>::single_empty_root() {
        if (single_empty_rooted_) {return root_;}

        // find all root nodes of the tree

        std::set<vertex_descriptor_t> old_roots;

        for (int i=0; i < int(tree_.size()); ++i) {
            old_roots.insert(i);
        }

        for (int i=0; i < int(tree_.size()); ++i) {
            for (auto &e: tree_.adj_edges(i)) {
                old_roots.erase( e.target() );
            }
        }

        if ( old_roots.size()==1 && tree_[ *old_roots.begin() ].cluster.empty() ) {
            root_ = *old_roots.begin();
        } else {
            // insert new root and point to all old roots
            root_ = tree_.add_vertex();

            for (auto old_root : old_roots) {
                tree_.add_edge(root_, old_root);
            }
        }

        single_empty_rooted_ = true;

        return root_;
    }

    template<class FunValue, class EvaluationPolicy>
    auto
    ClusterTree<FunValue,EvaluationPolicy>
    ::evaluate() {
        auto root = single_empty_root();

        if (!evaluated_) {
            dfs_evaluate(root);

            auto a = assignment_t(fn_.domains());
            evaluation_result_ = a.eval_determined(tree_[root].cluster.functions(), evaluation_policy_t());

            evaluated_ = true;
        }
        return evaluation_result_;
    }

    template<class FunValue, class EvaluationPolicy>
    auto
    ClusterTree<FunValue,EvaluationPolicy>::traceback() {

        assert(evaluated_);
        //assert(is_consistent());

        auto a = assignment_t(fn_.domains());

        dfs_traceback(single_empty_root(), a);

        return a;
    }

    template<class FunValue, class EvaluationPolicy>
    auto
    ClusterTree<FunValue,EvaluationPolicy>
    ::restricted_traceback(const std::set<var_idx_t> &variables,
        const assignment_t &assignment) {

        assert(evaluated_);

        //trace_variables = std::set<var_idx_t>(variables.begin(),variables.end());
        trace_variables = variables;

        auto a = assignment_t(assignment);

        dfs_traceback(single_empty_root(), a);

        trace_variables.reset();

        return a;
    }

    template<class FunValue, class EvaluationPolicy>
    auto
    ClusterTree<FunValue,EvaluationPolicy>::traceback_new(bool non_redundant, std::string non_redundant_mode) {

        assert(evaluated_);

        auto a = assignment_t(fn_.domains());

        if (non_redundant) {
            if (non_redundant_mode == "default") {
                
                if (!weights_tree_.defined()) {
                    // pass the necessary information from the cluster tree
                    // the first time the weights tree is used
                    std::vector<cluster_info_t> nodes;
                    vertex_descriptor_t p = 0;
                    std::vector<vertex_descriptor_t> preorder(tree_.size());
                    df_traversal(single_empty_root(), nodes, p, preorder);
                    weights_tree_.define(nodes, preorder);
                }

                dfs_traceback_new(single_empty_root(), a, true);
                weights_tree_.add_forbidden(a);
                return a;

            } else if (non_redundant_mode == "rejection") {
                dfs_traceback_new(single_empty_root(), a, false);
                
                for (const auto& fa : forbidden) {
                    bool sampled = true;
                    for (size_t i=0; i<a.size(); i++) {
                        if (a[i]!=fa[i]) {
                            sampled = false;
                            break;
                        }
                    }
                    if (sampled) {
                        throw std::exception();
                    }
                }

                forbidden.push_back(a);
                return a;

            }
        } else {
            dfs_traceback_new(single_empty_root(), a, false);
            return a;
        }
    }
}

#endif
