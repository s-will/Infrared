/** @brief Very simple python wrapper for treeedemposition using htd
 * Wrapper written by Sebastian Will, 2018
 *
 * This module supports computing a tree decomposition from a graph
 * INPUT: The input graph is specified by number of nodes and list of (hyper)-edges.
 * OUTPUT: The tree decomposition is obtained as list of bags and list of tree edges.
 *
 * ------------------------------
 * Example Python session:
 * ------------------------------
 * >>> import libhtdwrap as htd
 * >>> myhtd = htd.HTD(6,[[0,1,2,3],[4,5],[4,5,3,2],[1,5]])
 * >>> myhtd.decompose()
 * >>> myhtd.bags()
 * [[2, 3, 4, 5], [1, 2, 3, 5], [0, 1, 2, 3], [2, 3, 5], [1, 2, 3]]
 * >>> myhtd.edges()
 * [[0, 3], [3, 1], [1, 4], [4, 2]]
 */

#include <boost/python.hpp>
#include "boost_python_aux.hpp"

#include <htd/main.hpp>
#include <memory>
#include <stack>
#include <iterator>

namespace bpy = boost::python;

// helper to wrap plain pointer into unique ptr
template<class T>
auto wrap_unique_ptr(T *x) { return std::unique_ptr<T>(x); }

class HTD {
public:
    HTD(int num_vertices, const std::vector<std::vector<int>> &in_edges)
    {
        manager_ = wrap_unique_ptr(htd::createManagementInstance(htd::Id::FIRST));
        
        // Create a new graph instance which can handle (multi-)hyperedges.
        graph_ = wrap_unique_ptr(manager_->multiHypergraphFactory().createInstance());
        
        /*
         *  Add vertices to the sample graph.  The vertices of a graph are
         *  numbered in ascending order starting from 1.
         */
        graph_->addVertices(num_vertices);

        /* Add hyper edges */
        for (auto &edge : in_edges) {
            
            std::vector<htd::vertex_t> htd_edge;
            for (const auto x : edge ) {
                //make 1-based
                htd_edge.push_back(x+1);
            }
            
            graph_->addEdge( htd_edge );
        }
    };

    // vertices are index 0,...,num_vertices-1 in input and output
    void
    decompose() {
        // Get the default tree decomposition algorithm
        htd::ITreeDecompositionAlgorithm * algorithm =
            manager_->treeDecompositionAlgorithmFactory().createInstance();
        
        // call the tree decomposition algorithm
        htd::ITreeDecomposition * td =
            algorithm->computeDecomposition( *graph_.get() );

        // transform the td to lists of bags and edges
        // @todo: actually this conversion is a detour; consider to pass the td object directly and expose
        // methods root() and children()! Currently, the bags/edges representation requires to topo-sort
        // again before we can construct the cluster tree.
        //
        bags_.resize(td->vertexCount());
        
        auto stack = std::stack<htd::vertex_t>();
        stack.push(td->root());
        
        while(!stack.empty()) {
            auto v = stack.top();
            stack.pop();
            
             // make 0-based
            std::vector<int> bag;
            for(const auto x: td->bagContent(v)) { 
                bag.push_back(x-1);
            }
            bags_[v-1] = bag;

            const auto &cs = td->children(v);
            for (auto &c: cs) {
                //make 0-based (and cast to signed)
                edges_.push_back(std::vector<int>{ static_cast<int>(v-1), static_cast<int>(c-1) });
                stack.push(c);
            }
        }
    }

    //! @brief getter for the edges of the tree decomposition
    const auto &
    edges() {
        return edges_;
    }

    //! @brief getter for the bags of the tree decomposition
    const auto &
    bags() {
        return bags_;
    }

private:
    std::unique_ptr< htd::LibraryInstance > manager_;
    std::unique_ptr< htd::IMutableMultiHypergraph> graph_;

    using vec_vec_int_t = std::vector<std::vector<int>>;
    vec_vec_int_t edges_;
    vec_vec_int_t bags_;
};

BOOST_PYTHON_MODULE(libhtdwrap)
{
    register_vector_conversions<int>();
    register_vector_conversions<std::vector<int>>();

    bpy::class_<HTD, 
                boost::noncopyable>
        ("HTD", bpy::init<int, const std::vector<std::vector<int>>>())
        .def("decompose",&HTD::decompose)
        .def("edges",&HTD::edges,
             bpy::return_value_policy<bpy::copy_const_reference>())
        .def("bags",&HTD::bags,
             bpy::return_value_policy<bpy::copy_const_reference>())
        ;
}
