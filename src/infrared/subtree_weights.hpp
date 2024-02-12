#ifndef INFRARED_SUBTREE_WEIGHTS_HPP
#define INFRARED_SUBTREE_WEIGHTS_HPP

#include "graph.hpp"
#include "feature_network.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace ired {

    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy <FunValue>>
    class WeightsTree {

    public:
        using feature_network_t = FeatureNetwork<FunValue, EvaluationPolicy>;
        using assignment_t = typename feature_network_t::assignment_t;
        using cluster_t = typename feature_network_t::cluster_t;
        using evaluation_policy_t = typename feature_network_t::evaluation_policy_t;
        using fun_value_t = typename feature_network_t::fun_value_t;

        struct vertex_info_t {
            vertex_info_t()
                : Z(), w() {}
            vertex_info_t(fun_value_t Z_value, fun_value_t w_value)
                : Z(Z_value), w(w_value) {}
            fun_value_t Z;
            fun_value_t w; // TO DO: maek temporary
        };

        struct edge_info_t {
            edge_info_t(const assignment_t &assignment)
                : assignment(assignment) {}
            assignment_t assignment;
        };

        using tree_t = graph::adjacency_list<vertex_info_t, edge_info_t>;
        using vertex_descriptor_t = typename tree_t::vertex_descriptor_t;

        struct cluster_info_t {
            cluster_info_t(cluster_t cluster, std::vector<vertex_descriptor_t> children)
                : cluster(cluster), children(children) {}
            cluster_t cluster;
            std::vector<vertex_descriptor_t> children;
        };

        WeightsTree() {}

        ~WeightsTree() {}

        bool
        defined() {
            return defined_;
        }

        void
        define(std::vector<cluster_info_t> &nodes, std::vector<vertex_descriptor_t> &preorder) {
            vertex_descriptor_t root_ = tree_.add_vertex();
            path_.push_back(root_);
            preorder_ = preorder;
            cluster_nodes_ = nodes;
            defined_ = true;
        }

        void
        add_forbidden(assignment_t &a) {
            update_tree_structure(a);
            update_weights(a);
            path_.erase(path_.begin() + 1, path_.end());

            std::string filename = "weights_graph_" + std::to_string(count) + ".dot";
            record_tree(filename);
        }

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

        std::vector<vertex_descriptor_t> path_;
        std::vector<cluster_info_t> cluster_nodes_;
        std::vector<vertex_descriptor_t> preorder_;

        bool defined_ = false;

        int count = 0;

        void
        add_node(const edge_info_t& e) {
            vertex_descriptor_t current_node = path_.back();
            vertex_descriptor_t new_node = tree_.add_vertex();

            tree_[new_node].Z = evaluation_policy_t::zero();
            tree_[new_node].w = evaluation_policy_t::zero();
            
            tree_.add_edge(current_node, new_node, e);
            path_.push_back(new_node);
        }

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