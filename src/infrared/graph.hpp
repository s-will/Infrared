#ifndef GRAPH_HPP
#define GRAPH_HPP

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
 * @brief graph data structure(s)
 */

namespace graph {
    /** @brief Graph in adjacency list representation
     */
    template<class Vertex, class Edge>
    class adjacency_list {
    public:
        using vertex_descriptor_t = size_t;

        using vertex_t = Vertex;
        using edge_info_t = Edge;

        class edge_t: public edge_info_t {
        public:
            edge_t(vertex_descriptor_t source,
                   vertex_descriptor_t target,
                   const edge_info_t &edge_info)
                : source_(source),
                  target_(target),
                  edge_info_t(edge_info) {
            }

            auto
            source() const {
                return this->source_;
            }

            auto
            target() const {
                return this->target_;
            }

        private:
            vertex_descriptor_t source_;
            vertex_descriptor_t target_;
        };

        //!@brief constructor
        adjacency_list() {
        }

        //!@brief add vertex
        vertex_descriptor_t
        add_vertex(const vertex_t &v) {
            auto idx = vertices_.size();
            vertices_.push_back(v);
            adj_edges_.push_back(std::vector<edge_t>());
            return idx;
        }

        //!@brief add vertex
        vertex_descriptor_t
        add_vertex() {
            return add_vertex(vertex_t());
        }

        //!@brief add egde
        auto &
        add_edge(vertex_descriptor_t source,
                 vertex_descriptor_t target,
                 const edge_info_t &e) {
            assert(source < vertices_.size());
            adj_edges_[source].push_back(edge_t(source,target,e));
            return adj_edges_[source].back();
        }

        //!@brief add egde
        auto &
        add_edge(vertex_descriptor_t source,
                 vertex_descriptor_t target) {
            return add_edge(source, target, edge_info_t());
        }

        const auto &
        operator [] (vertex_descriptor_t i) const {
            return vertices_[i];
        }

        auto &
        operator [] (vertex_descriptor_t i) {
            return vertices_[i];
        }

        //!@brief iterable of (descriptors of) adjacent edges
        const auto &
        adj_edges(vertex_descriptor_t i) const {
            return adj_edges_[i];
        }

        //!@brief iterable of (descriptors of) adjacent edges
        auto &
        adj_edges(vertex_descriptor_t i) {
            return adj_edges_[i];
        }

        //!@brief number of vertices
        auto
        size() {
            return vertices_.size();
        }

    private:
        //! vector of vertices
        std::vector<vertex_t> vertices_;

        //! adjacency lists
        std::vector<std::vector<edge_t>> adj_edges_;
    }; // end class adjacency_list

    // Example code for dfs; implementing visitors (like e.g. in boost::graph)
    // seems to be overkill for our simple purposes
    //
    // template<class Graph>
    // void
    // dfs(const Graph &g, Graph::vertex_descriptor_t root, std::set<Graph::vertex_descriptor_t> &visited) {
    //     if(visited.in(root)) {
    //         return;
    //     }
    //     visited.add(root);
    //
    //     for(auto e: g.adj_edges(root)) {
    //         examine_edge(g,e);
    //         dfs(g, e.target(), visited);
    //         finish_edge(g,e);
    //     }
    // };
} //end namespace graph
#endif // GRAPH_HPP
