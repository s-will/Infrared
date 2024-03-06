#ifndef INFRARED_SUBTREE_WEIGHTS_HPP
#define INFRARED_SUBTREE_WEIGHTS_HPP

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
 * @brief Defines the weights tree.
 */

#include "graph.hpp"
#include "feature_network.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace ired {
    /**
     * @brief Tree data structure containing information about weights 
     * that should be added based on previously generated assignments 
     * in order to achieve non-redundancy
     */
    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy <FunValue>>
    class WeightsTree {

    public:
        using feature_network_t = FeatureNetwork<FunValue, EvaluationPolicy>;
        using assignment_t = typename feature_network_t::assignment_t;
        using cluster_t = typename feature_network_t::cluster_t;
        using evaluation_policy_t = typename feature_network_t::evaluation_policy_t;
        using fun_value_t = typename feature_network_t::fun_value_t;
        using function_t = typename feature_network_t::function_t;
        using var_idx_t = typename feature_network_t::var_idx_t;
        using var_value_t = int;

        //! @brief information at a vertex: weight Z
        struct vertex_info_t {
            vertex_info_t(fun_value_t Z_value)
                : Z(Z_value) {}
            fun_value_t Z;
        };

        //! @brief information at an edge: partial assignment a
        struct edge_info_t {
            edge_info_t(std::vector<var_value_t> &assignment)
                : assignment(assignment) {}
            std::vector<var_value_t> assignment;
        };

        using tree_t = graph::adjacency_list<vertex_info_t, edge_info_t>;
        using vertex_descriptor_t = typename tree_t::vertex_descriptor_t;

        // information about a node from the cluster tree:
        // cluster (variables, constraints and functions) + list of children
        struct cluster_info_t {
            cluster_info_t(std::vector<const function_t *> functions, size_t num_children)
                : functions(functions), num_children(num_children) {}
            std::vector<const function_t*> functions;
            size_t num_children;
        };

        WeightsTree() {}

        ~WeightsTree() {}

        //! @brief checks whether the necessary information from the cluster tree has been passed
        bool
        defined() {
            return defined_;
        }

        //! @brief sets the cluster tree nodes information and their order in df traversal
        void
        define(std::vector<std::vector<var_idx_t>> &diffs, 
               std::vector<cluster_info_t> &nodes, 
               std::vector<vertex_descriptor_t> &preorder) {
            diffs_ = diffs;
            cluster_nodes_ = nodes;
            preorder_ = preorder;

            vertex_descriptor_t root_ = tree_.add_vertex(vertex_info_t(evaluation_policy_t::zero()));
            path_.push_back(root_);
            w_.push_back(evaluation_policy_t::one());

            defined_ = true;
        }

        void
        clear() {
            path_.erase(path_.begin() + 1, path_.end());
            w_.erase(w_.begin() + 1, w_.end());
            w_[0] = evaluation_policy_t::one();
        }

        /**
         * @brief Add the necessary nodes and edges to represent a new assignment
         *
         * @param a partial assignment
         * @param diff a list of newly assigned variables
         *
         * At each level of the weights tree the new assignment is compared to 
         * the existing partial assignments. If there is a match, the search 
         * continues in the respective subtree. Otherwise a new node is added.
         */
        void
        update_tree_structure(assignment_t &a, std::vector<var_idx_t> &diff) {
            bool exists = false;
            for (auto &e: tree_.adj_edges(path_.back())) {
                bool match = true;
                for (size_t i=0; i<diff.size(); ++i) {
                    if (a[diff[i]] != e.assignment[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    path_.push_back(e.target());
                    w_.push_back(funs_product(cluster_nodes_[path_.size()-1], a));
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                add_node(a, diff);
            }
        }

        /**
         * @brief Update the weights after adding a new assignment
         *
         * @param a total assignment
         *
         * Following the assignment path from the leaf to the root, 
         * update each weight by adding the product of all functions in the subtree 
         */
        void
        update_weights(assignment_t &a) {
            tree_[path_.back()].Z = w_.back();
            for (size_t i = path_.size()-2; i != SIZE_MAX; --i) {
                w_[i] = evaluation_policy_t::mul(w_[i], w_[i+1]);
                tree_[path_[i]].Z = evaluation_policy_t::plus(tree_[path_[i]].Z, w_[i]);
            }
            //record_tree();
        }

        /**
         * @brief Get specific subtree weight
         *
         * @param v cluster tree node
         * @param a assignment
         */
        fun_value_t
        get_weight(vertex_descriptor_t v, const assignment_t &a) {
            vertex_descriptor_t current_node = path_.back();
            for (size_t k=path_.size()-1; k<preorder_[v]; ++k) {
                auto diff = diffs_[k];
                bool exists = false;
                for (auto &e: tree_.adj_edges(current_node)) {
                    bool match = true;
                    for (size_t i=0; i<diff.size(); ++i) {
                        if (a[diff[i]] != e.assignment[i]) {
                            match = false;
                            break;
                        }
                    }
                    if (match) {
                        current_node = e.target();
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    return evaluation_policy_t::zero();
                }
            }
            return tree_[current_node].Z;
        }

        fun_value_t
        get_value() {
            return w_.back();
        }

        fun_value_t
        get_accumulated_weight() {
            return tree_[0].Z;
        }

    private:
        tree_t tree_;

        std::vector<vertex_descriptor_t> path_;
        std::vector<fun_value_t> w_;
        
        std::vector<cluster_info_t> cluster_nodes_;
        std::vector<vertex_descriptor_t> preorder_;
        std::vector<std::vector<var_idx_t>> diffs_;

        bool defined_ = false;

        /**
         * @brief Add new node to the weights tree
         *
         * @param e assignment for the connecting edge
         *
         * Every new node is initialized with weight 0.
         */
        void
        add_node(assignment_t &a, std::vector<var_idx_t> &diff) {
            vertex_descriptor_t current_node = path_.back();
            
            vertex_descriptor_t new_node = tree_.add_vertex(vertex_info_t(evaluation_policy_t::zero()));
            path_.push_back(new_node);
            
            fun_value_t product = funs_product(cluster_nodes_[path_.size()-1], a);
            w_.push_back(product);

            std::vector<var_value_t> e;
            for (size_t i=0; i<diff.size(); ++i) {
                e.push_back(a[diff[i]]);
            }
            tree_.add_edge(current_node, new_node, e);
        }

        /**
         * @brief Compute the product of cluster functions for a given assignment
         *
         * @param cluster_node cluster information (cluster + list of children)
         * @param a total assignment
         * 
         * @return value of the product
         */
        auto
        funs_product(cluster_info_t &cluster_node, assignment_t &a) {
            auto funs = cluster_node.functions;
            size_t num_funs = funs.size() - cluster_node.num_children;

            auto x = evaluation_policy_t::one();
            for (size_t i = 0; i < num_funs; ++i) {
                x = evaluation_policy_t::mul(x, (*funs[i])(a));
            }
            return x;
        }

        /**
         * @brief Store the current state of the weights tree in a dotfile
         */
        void
        record_tree() {
            std::ofstream dotFile("weights_graph.dot");
            dotFile << "graph G {" << std::endl;
            for (vertex_descriptor_t v = 0; v < tree_.size(); ++v) {
                dotFile << "  var" << v << " [label=\"" << tree_[v].Z << "\"];" << std::endl;
            }
            for (vertex_descriptor_t v = 0; v < tree_.size(); ++v) {
                for(const auto &e: tree_.adj_edges(v)) {
                    std::string a = "";
                    for (size_t i=0; i<e.assignment.size(); ++i) {
                        a += std::to_string(e.assignment[i]);
                    }
                    dotFile << "  var" << e.source() << " -- var" << e.target() << " [label=\"" << a << "\"];" << std::endl;           
                }
            }
            dotFile << "}" << std::endl;
            dotFile.close();
        }
    };

}

#endif