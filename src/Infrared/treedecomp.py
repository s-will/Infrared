#!/usr/bin/env python3

#
# InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

## @file
##
## Storing and manipulating tree decompositions
##
## Supports computation by delegate
##
## * generate tree decompositions usng libhtd (via python wrapper) or calling TDlib
##
## * plot dependency graphs and tree decompositions using dot
##
## NOTE: TDlib writes files to disk, currently it is *not* safe to use
## in distributed computations. In contrast, libhtd (default) does not
## write to disk at all
##
## For using TDlib, it has to be in the Java class path (environment
## variable CLASSPATH). The default tree decomposition by libhtd needs
## to find the shared library libhtd.so; this could require to set
## LD_LIBRARY_PATH manually.

import os
import subprocess
import re
from math import sqrt,ceil
import itertools
import abc
import random

def seed(seed=None):
    """
    @brief seed treedecomposition random number generator

    Seeds the python random number generator (globally!).
    Without argument or seed==None, use pythons built-in random.seed() to generate
    a seed.
    """
    random.seed(seed)

## @brief Class to hold a tree decomposition
##
## contains
##
## * bags   a list of 0-based indexed lists of vertices of the decomposed graph (or "variables")
##
## * edges  a list of edges; an edge is a list of bags (bag1, bag2)
##
## * adj    adjacency lists (as defined by edges)
##
## edges are directed, an edge from bag x to y is represented by (x,y)
class TreeDecomposition:

    ## @brief construct from bags and edges
    ##
    ## @param bags lists of indices of the variables (vertices) in the
    ## bag
    ##
    ## @param edges list of edges between the bags; edges are
    ## specified as pairs (x,y), where x is the index of the parent
    ## bag and y, of the child bag (0-based indices according to
    ## bags).
    ##
    ## Edges must be directed from root(s) to leaves. Typically the
    ## bags and edges result from calling a tree decomposition
    ## algorithm. The supported algorithms return correctly oriented
    ## edges.
    def __init__(self, bags, edges):
        # copy bags
        self._bags = list(bags)
        # ensure that all edges are tuples or lists of length 2
        assert(all(len(x)==2 for x in edges))
        # copy edges and convert all edges to pairs (even if they have been lists)
        self._edges = list(map(lambda x: (x[0],x[1]), edges))

        self.update_adjacency_lists()

    ## @brief list of bags
    def get_bags(self):
        return self._bags

    ## @brief list of edges
    def get_edges(self):
        return self._edges

    ## @brief Comute adjacency list representation
    ##
    ## @param n number of nodes
    ## @param edges list of edges
    @staticmethod
    def adjacency_lists(n, edges):
        adj = { i:[] for i in range(n) }
        for (i,j) in edges:
            adj[i].append(j)
        return adj

    ## @brief Update the adjacency representation
    ##
    ## Updates the adjacency representation in adj according to the
    ## number of bags and edges
    def update_adjacency_lists(self):
        self.adj = self.adjacency_lists(len(self._bags),self._edges)

    ## @brief Toppological sort of bags
    ##
    ## @returns sorted list of bag indices
    def toposorted_bag_indices(self):
        n = len(self._bags)

        visited = set()
        sorted = list()

        def toposort_component(i):
            visited.add(i)
            for j in self.adj[i]:
                if not j in visited:
                    toposort_component(j)
            sorted.append(i)

        for i in range(n):
            if not i in visited:
                toposort_component(i)
        return sorted[::-1]

    ## @brief Difference set
    ##
    ## @param xs first list
    ## @param ys second list
    ##
    ## @return ys setminus xs
    ##
    ## For bags xs and ys, this computes the introduced variable
    ## indices when going from parent xs to child ys.
    @staticmethod
    def diff_set(xs,ys):
        return [ y for y in ys if y not in xs ]

    ## @brief Separator set
    ##
    ## @param xs first list
    ## @param ys second list
    ##
    ## @return overlap of xs and ys
    ##
    ## For bags xs and ys, this computes the 'separator' of xs and ys,
    ## i.e. the common variables.
    @staticmethod
    def sep_set(xs,ys):
        return [ y for y in ys if y in xs ]

    ## @brief Get tree width
    def treewidth(self):
        return max([len(bag) for bag in self._bags]) - 1

    ## @brief Write tree decomposition in dot format
    ## @param out output file handle
    def writeTD(self, out):
        def baglabel(bag):
            if len(bag)==0:
                return ""

            lwidth = ceil( sqrt(len(bag)) )
            lnum   = ceil( len(bag) / lwidth )
            xs = [str(i) for i in sorted(bag)]
            lines=list()
            for i in range(0,lnum):
                lines.append(" ".join(xs[i*lwidth:(i+1)*lwidth]))
            return "\\n".join(lines)

        out.write("digraph G {\n\n")

        for bagid,bag in enumerate(self._bags):
            label = baglabel(bag)
            out.write( "\tbag{} [label=\"{}\"]\n".format(bagid+1,label) )

        out.write("\n\n")

        for (x,y) in self._edges:
            edgelabel = " ".join( [ str(x) for x in self.diff_set(self._bags[x],self._bags[y] )] )
            out.write( "\tbag{} -> bag{}  [label=\"{}\"]\n".format(x+1,y+1,edgelabel) )

        out.write("\n}\n")

    ## @brief Guarantee a certain maximal diff set size
    ##
    ## @param maxdiffsize maximum size of diff sets after transformation
    ##
    ## @see diff_set()
    ##
    ## @pre the tree is connected
    ##
    ## Expands the tree by inserting a minimum number of in-between
    ## bags whereever the diff set size exceeds maxdiffsize. This can
    ## limit the complexity of certain algorithm based on the tree (in
    ## particular sampling in Infrared).
    ##
    def expand_treedecomposition(self, maxdiffsize=1):

        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in range(0, len(l), n):
                yield l[i:i + n]

        n = len(self._bags)
        root = self.toposorted_bag_indices()[0]

        next_bag_idx = n

        ## perform depth first traversal
        visited = set()
        stack = [(None,root)] # parent idx, bag index
        while stack:
            #pop
            (u,v) = stack[-1]
            stack = stack[:-1]

            if v in visited: continue
            visited.add(v)

            # push children on stack
            for w in self.adj[v]:
                stack.append((v,w))

            if u is not None:
                # determine diff set
                diff = self.diff_set(self._bags[u],self._bags[v])
                if len(diff) > maxdiffsize:
                    sep = self.sep_set(self._bags[u],self._bags[v])

                    if (u,v) in self._edges:
                        self._edges.remove((u,v))

                    if (v,u) in self._edges:
                        self._edges.remove((v,u))

                    last_bag_idx = u

                    newbag = sep

                    for ext in chunks(diff[:-1],maxdiffsize):
                        newbag.extend(ext)
                        self._bags.append(newbag[:])
                        self._edges.append([ last_bag_idx, next_bag_idx])
                        last_bag_idx = next_bag_idx
                        next_bag_idx += 1

                    self._edges.append([ last_bag_idx, v])

        self.update_adjacency_lists()

# ##########################################################
# writing graphs to files in different formats
#

## @brief Plot graph
##
## @param graphfile file of graph in dot format
##
## The graph is plotted and written to a pdf file by calling
## graphviz's dot tool on th dot file.
def dotfile_to_pdf(graphfile):
    outfile = re.sub(r".dot$",".pdf",graphfile)
    subprocess.check_output(["dot","-Tpdf","-o",outfile,graphfile])

def dotfile_to_png(graphfile):
    outfile = re.sub(r".dot$",".png",graphfile)
    subprocess.check_output(["dot","-Tpng","-o",outfile,graphfile])

## @brief Write graph in dgf format
##
##@param out output file handle
##@param num_nodes number of nodes
##@param edges the edges of the graph
def write_dgf(out, num_nodes, edges):
    edge_num=len(edges)

    out.write("p tw {} {}\n".format(num_nodes, edge_num))
    for (u,v) in sorted(edges):
        out.write("e {} {}\n".format(u+1,v+1))

## @brief Write graph in dot format
##
## @param out output file handle
## @param num_nodes number of nodes
## @param edges the edges of the graph
def write_dot(out, num_nodes, edges):
    edge_num=len(edges)

    out.write("graph G{\n\n")

    for v in range(num_nodes):
        out.write("\tnode{idx} [label=\"{idx}\"]\n".format(idx=v))

    for (u,v) in edges:
        out.write("\tnode{} -- node{}\n".format(u,v))

    out.write("\n}\n")


# ##########################################################
# Interface TDlib
#
# specific functions to call TDlib's tree decomposition tool
# and parse the result
#

## @brief Parse tree decomposition as written by TDlib
##
## @param tdfh file handle to tree decomposition in dot format
## @param num_nodes number of nodes
## @returns tree decomposition
##
## Assume file is in td format as written by the tool TDlib.
def parseTD_TDlib(tdfh, num_nodes):

    bags = list()
    edges = list()

    for line in tdfh:
        # catch bags "New_vertex"
        if re.search("New_Vertex",line):
            bags.append([])
            continue

        m = re.search(r"bag(\d+).+label=\"(.+)\"",line)
        if m:
            bagid=int(m.group(1))
            label = m.group(2)
            labels = re.findall('\d+', label)
            labels = [ int(label) for label in labels ]
            if bagid!=len(bags)+1:
                raise IOError("Bag indices in td file must be consecutive (at bag {})!".format(bagid))
            bags.append(labels)
        else:
            m = re.search(r"bag(\d+) -- bag(\d+)", line)
            if m:
                edges.append((int(m.group(1)),
                              int(m.group(2))))
    # decrease bag labels

    def dec(xs):
        return [ [x-1 for x in ys]  for ys in xs]

    bags = dec(bags)
    edges = dec(edges)

    # add missing bags
    present = set()
    for b in bags:
        for x in b:
            present.add(x)

    for i in range(num_nodes):
        if not i in present:
            bags.append([i])

    return TreeDecomposition(bags,edges)

## @brief Compute tree decomposition of a graph by TDlib
##
## Generates tree decomposition and writes the result to file.
##
## @param filename base name for input .dfg and output .td files
## @param num_nodes number of nodes
## @param edges specifies edges of a graph; nodes are indexed 0-based
## @param strategy integer code of decomposition strategy:
##   1) Permutation To Tree Decomposition with GreedyDegree
##   2) Permutation To Tree Decomposition with GreedyFillIn
##   3) Permutation To Tree Decomposition with Lex-BFS: Lexicographic Breadth First Search
##   4) Permutation To Tree Decomposition with MCS; Maximum Cardinality Search Algorithm
##   5) PermutationGuesser
def makeTDFile_TDlib( filename, num_nodes, edges,
                      *, strategy=2 ):

    inname = filename+".dgf"
    outname = filename+".td"

    dotfilename = filename+".dot"

    with open(inname,"w") as dgf:
        write_dgf(dgf,num_nodes,edges)
        dgf.close()

        cmd = ["java", # "-cp", tdlib_home,
               "nl.uu.cs.treewidth.TreeDecomposer",
               str(strategy), inname, outname
               #, filename+".dot"
        ]
        try:
            subprocess.check_output(cmd)
        except:
            print("ERROR: Could not call TDlib; please make sure that TDlib is installed in the Java class path")
            print("       Possibly, you just need to set the environment variable CLASSPATH, like")
            print("       export CLASSPATH=$tdlib_home/treewidth-java")
            exit(-1)

        os.remove(inname)

        return outname

# ###########################################################
# Interface tree demposition libraries
#

## @brief Base class of tree decomposition factories
##
## A TD factory needs to provide a method create to produce a class TreeDecomposition
## given the number of variables and the list of dependencies;
## dependencies are lists of lists of 0-based indices of the
## variables that respectively depend on each other
##
class TreeDecompositionFactoryBase:
    def __init__(self):
        pass

    ## @brief Create tree decomposition
    ##
    ## @param size number of nodes in the dependency graph
    ## @param dependencies specifies edges of the dependency (hyper-)graph; nodes are indexed 0-based
    ##
    ## @return tree decomposition (object of TreeDecomp)
    @abc.abstractmethod
    def create(self, size, dependencies):
        return

    ## @brief Expand non-binary dependencies to cliques of binary deps
    ## @param dependencies list of dependencies
    ## @return list of binary dependencies
    @staticmethod
    def expand_to_cliques(dependencies):
        bindeps = list()
        for d in dependencies:
            bindeps.extend( itertools.combinations(d,2) )
        return bindeps


## @brief Tree decomposition factory using HTD
class HTDTreeDecompositionFactory(TreeDecompositionFactoryBase):
    def __init__(self, maxdiffsize=1):
        self.maxdiffsize = maxdiffsize
        pass

    ## @brief Create tree decomposition
    def create(self, size, dependencies):
        bindependencies  = self.expand_to_cliques(dependencies)

        from libhtdwrap import HTD
        myhtd = HTD(size,bindependencies)
        myhtd.decompose()

        td = TreeDecomposition(myhtd.bags(), myhtd.edges())
        td.expand_treedecomposition(self.maxdiffsize)

        return td

## @brief Tree decomposition factory using TDlib
class TDLibTreeDecompositionFactory(TreeDecompositionFactoryBase):
    ## @brief construct
    ## @param strategy TDlib's strategy (see help page of TDlib)
    ## @param tmpfile file for tdlib's output
    ## @todo automatically generate unique / thread-safe tmp name
    def __init__(self, strategy=2, tmpfile="/tmp/treedecomp-tmp~"):
        self.strategy = strategy
        pass

    ## @brief Make tree decomposition using TDlib
    ##
    ##
    ## @return tree decomposition (object of TreeDecomp)
    ##
    def create(self, size, dependencies):
        bindependencies = self.expand_to_cliques(dependencies)

        makeTDFile_TDlib( self.tmpfile, size, bindependencies, strategy=self.strategy )
        tdfile = self.tmpfile+".td"
        with open(tdfile) as tdfh:
            td =  parseTD_TDlib( tdfh, size )
        os.remove(tdfile)
        return td

## @brief Tree decomposition factory using networkx
class NXTreeDecompositionFactory(TreeDecompositionFactoryBase):
    def __init__(self, maxdiffsize=1):
        self.maxdiffsize = maxdiffsize
        self.iterations = 20
        pass

    ## @brief Create tree decomposition
    def create(self, size, dependencies):

        # iterate tree decomposition heuristics over randomized graphs

        def translate( xs, perm ):
            return [ [ perm[y] for y in ys ] for ys in xs ]

        bindependencies = self.expand_to_cliques(dependencies)

        # produce networkx graph from size, dependencies
        from networkx import Graph

        from networkx.algorithms.approximation.treewidth import treewidth_min_fill_in, treewidth_min_degree

        import random

        perm = list( range(size) )

        best_width = None

        for k in range(self.iterations):
            random.shuffle( perm )
            inv_perm = { perm[i]:i for i in perm }

            G = Graph()
            G.add_nodes_from( range(size) )
            G.add_edges_from( translate( bindependencies, perm ) )

            width, tree = treewidth_min_fill_in(G)

            if best_width is None or width < best_width:
                 best_width, best_tree = width, tree
                 best_inv_perm = inv_perm

        bags = list(map(list, best_tree.nodes))
        edges = [(bags.index(list(i)),bags.index(list(j))) for i,j in best_tree.edges]
        bags = translate( bags, best_inv_perm )
        td = TreeDecomposition(bags, edges)
        td.expand_treedecomposition(self.maxdiffsize)

        return td

## @brief default tree decomposition factory
TreeDecompositionFactory = NXTreeDecompositionFactory

## @brief the available, predefined td factories with some description
_td_factories = [
        ["nx", "using networkx module", NXTreeDecompositionFactory],
        ["tdlib", "using TDlib, strategy 2", lambda: TDLibTreeDecompositionFactory(strategy=2)],
        ["htd", "using libhtd", HTDTreeDecompositionFactory]
        ]

## @brief get a tree decomposition factory by a descriptor string
def td_factory_from_descriptor(descriptor):
    for x in _td_factories:
        if descriptor == x[0]:
            return x[2]()
    return None

def get_td_factory_descriptors():
    return [ x[0] for x in _td_factories ]

# End of Interface tree demposition libraries
# ----------------------------------------------------------

if __name__ == "__main__":
    pass
