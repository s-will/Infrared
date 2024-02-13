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

        //! @brief information at a vertex: weight Z 
        struct vertex_info_t {
            vertex_info_t()
                : Z(), w() {}
            vertex_info_t(fun_value_t Z_value, fun_value_t w_value)
                : Z(Z_value), w(w_value) {}
            fun_value_t Z;
            fun_value_t w;  //TO DO: make temporary
        };

        //! @brief information at an edge: partial assignment a
        struct edge_info_t {
            edge_info_t(const assignment_t &assignment)
                : assignment(assignment) {}
            assignment_t assignment;
        };

        using tree_t = graph::adjacency_list<vertex_info_t, edge_info_t>;
        using vertex_descriptor_t = typename tree_t::vertex_descriptor_t;

        // information about a node from the cluster tree:
        // cluster (variables, constraints and functions) + list of children
        struct cluster_info_t {
            cluster_info_t(cluster_t cluster, std::vector<vertex_descriptor_t> children)
                : cluster(cluster), children(children) {}
            cluster_t cluster;
            std::vector<vertex_descriptor_t> children;
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
        define(std::vector<cluster_info_t> &nodes, std::vector<vertex_descriptor_t> &preorder) {
            vertex_descriptor_t root_ = tree_.add_vertex();
            path_.push_back(root_);
            preorder_ = preorder;
            cluster_nodes_ = nodes;
            defined_ = true;
        }

        /**
         * @brief Modify the tree with a new assignment
         *
         * @param a generated assignment by traceback
         */
        void
        add_forbidden(assignment_t &a) {
            update_tree_structure(a);
            update_weights(a);
            path_.erase(path_.begin() + 1, path_.end());

            std::string filename = "weights_graph_" + std::to_string(count) + ".dot";
            record_tree(filename);
        }

        /**
         * @brief Get specific subtree weight
         *
         * @param v cluster tree node
         * @param a assignment
         */
        fun_value_t
        get_weight(vertex_descriptor_t v, assignment_t &a) {
            vertex_descriptor_t current_node = path_[0];
            for (size_t i = 0; i < preorder_[v]; ++i) {
                auto& parent = (cluster_nodes_[i]).cluster;
                auto& child = cluster_nodes_[i+1].cluster;

                auto diff = child.diff_vars(parent);

                bool exists = false;
                for (auto &e: tree_.adj_edges(current_node))  {
                    if (std::all_of(std::begin(diff), std::end(diff), [&](const auto &i) {
                            return a[i] == e.assignment[i];
                        })) {
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

    private:
        feature_network_t fn_;
        tree_t tree_;

        // stores the tree nodes from the weights tree 
        // corresponding to the currently proccessed assignment
        std::vector<vertex_descriptor_t> path_;
        
        std::vector<cluster_info_t> cluster_nodes_;
        std::vector<vertex_descriptor_t> preorder_;

        bool defined_ = false;

        int count = 0; //used for generating a diferent dotfile name

        /**
         * @brief Add new node to the weights tree
         *
         * @param e assignment for the connecting edge
         *
         * Every new node is initialized with weight 0.
         */
        void
        add_node(const edge_info_t& e) {
            vertex_descriptor_t current_node = path_.back();
            vertex_descriptor_t new_node = tree_.add_vertex();

            tree_[new_node].Z = evaluation_policy_t::zero();
            tree_[new_node].w = evaluation_policy_t::zero();
            
            tree_.add_edge(current_node, new_node, e);
            path_.push_back(new_node);
        }

        /**
         * @brief Add the necessary nodes and edges to represent a new assignment
         *
         * @param a total assignment
         *
         * At each level of the weights tree the new assignment is compared to 
         * the existing partial assignments. If there is a match, the search 
         * continues in the respective subtree. Otherwise a new node is added.
         */
        void
        update_tree_structure(assignment_t &a) {
            for (size_t i = 0; i < cluster_nodes_.size()-1; ++i) {
                auto& parent = (cluster_nodes_[i]).cluster;
                auto& child = cluster_nodes_[i+1].cluster;

                auto diff = child.diff_vars(parent);

                bool exists = false;
                for (auto &e: tree_.adj_edges(path_.back()))  {
                    if (std::all_of(std::begin(diff), std::end(diff), [&](const auto &i) {
                            return a[i] == e.assignment[i];
                        })) {
                        path_.push_back(e.target());
                        exists = true;
                        break;
                    }
                }

                if (!exists) {
                    add_node(a);
                }
            }
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
            auto funs = cluster_node.cluster.functions();
            size_t num_funs = funs.size() - cluster_node.children.size();

            auto x = evaluation_policy_t::one();
            for (size_t i = 0; i < num_funs; ++i) {
                x = evaluation_policy_t::mul(x, (*funs[i])(a));
            }
            return x;
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
            tree_[path_.back()].w = funs_product(cluster_nodes_.back(), a);
            tree_[path_.back()].Z = tree_[path_.back()].w;
            
            for (size_t i = path_.size() - 2; i > 0; --i) {
                auto w = funs_product(cluster_nodes_[i], a);
                for(const auto &child: cluster_nodes_[i].children) {
                    w = evaluation_policy_t::mul(w, tree_[path_[preorder_[child]]].w);
                }

                tree_[path_[i]].w = w;
                tree_[path_[i]].Z = evaluation_policy_t::plus(tree_[path_[i]].Z, w);
            }
        }

        /**
         * @brief Store the current state of the weights tree in a dotfile
         *
         * @param filename name of the dotfile
         */
        void
        record_tree(std::string filename) {
            std::ofstream dotFile(filename);
            dotFile << "graph G {" << std::endl;
            for (vertex_descriptor_t v = 0; v < tree_.size(); ++v) {
                dotFile << "  var" << v << " [label=\"" << tree_[v].Z << ", " << tree_[v].w << "\"];" << std::endl;
            }
            for (vertex_descriptor_t v = 0; v < tree_.size(); ++v) {
                for(const auto &e: tree_.adj_edges(v)) {
                    std::string a = "";
                    for (size_t i=0; i<e.assignment.size(); ++i) {
                        a += "x" + std::to_string(i) + "=" + std::to_string(e.assignment[i]) + ", ";
                    }
                    dotFile << "  var" << e.source() << " -- var" << e.target() << " [label=\"" << a << "\"];" << std::endl;           
                }
            }
            dotFile << "}" << std::endl;
            dotFile.close();
            count ++;
        }
    };

}

#endif