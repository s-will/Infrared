#!/usr/bin/env python3

##############################
# Compute (by delegate) and handle tree decompositions
# * generate tree decompositions usng libhtd (via python wrapper) or calling TDlib
# * plot dependency graphs and tree decompositions using dot
#
# NOTE: in constrast to TDlib, libhtd does not need to write intermediary files to disk
# ATTENTION: TDlib is *not* safe to use in distributed computations (@todo)
#
# for using TDlib, the Java class path needs to include the library!

import os
import subprocess
import re
from math import sqrt,ceil

import libhtdwrap as htd

############################################################
##
## Storing and manipulating tree decompositions
##

## Class to hold a tree decomposition
##
## contains
## bags   a list of 0-based indexed lists of vertices of the decomposed graph (or "variables")
## edges  a list of edges; an edge is a list of bags (bag1, bag2)
## adj    adjacency lists (as defined by edges)
##
## edges are directed, an edge from bag x to y is represented by (x,y)
##
class TreeDecomp:
    def __init__(self,bags,edges):
        # copy bags
        self.bags = list(bags)
        # ensure that all edges are tuples or lists of length 2
        assert(all(len(x)==2 for x in edges))
        # copy edges and convert all edges to pairs (even if they have been lists)
        self.edges = list(map(lambda x: (x[0],x[1]), edges))

        self.update_adjacency_lists()

    @staticmethod
    def adjacency_lists(n,edges):
        """
        compute adjacency lists
        """
        adj = { i:[] for i in range(n)}
        for (i,j) in edges:
            adj[i].append(j)
        return adj

    def update_adjacency_lists(self):
        self.adj = self.adjacency_lists(len(self.bags),self.edges)

    def toposorted_bag_indices(self):
        """
        @brief perform topological sort
        @param n number of nodes
        @param adj adjacency lists
        @returns sorted nodes
        """
        n = len(self.bags)

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

    # difference set from bag xs to bag ys
    @staticmethod
    def diff_set(xs,ys):
        """
        Compute directed difference between bags ys-xs
        """
        return [ y for y in ys if y not in xs ]

    # seperator set between bags xs and ys
    @staticmethod
    def sep_set(xs,ys):
        return [ y for y in ys if y in xs ]

    @staticmethod
    def max_node(edges):
        """
        Maximum node occuring in a list of edges
        """
        maxnode=0
        if len(edges)>0:
            maxnode=max([ max(u,v) for (u,v) in edges ])
            return maxnode

    def treewidth(self):
        """
        Get tree width
        """
        return max([len(bag) for bag in self.bags]) - 1

    # def inc_edges():
    #     """Make 0-based edges 1-based"""
    #     self.edges = [ (i+1,j+1) for (i,j) in self.edges ]
    #

    def writeTD(self, out):
        """
        Write tree decomposition in dot format
        @param out output file handle
        """

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

        for bagid,bag in enumerate(self.bags):
            label = baglabel(bag)
            out.write( "\tbag{} [label=\"{}\"]\n".format(bagid+1,label) )

        out.write("\n\n")

        for (x,y) in self.edges:
            edgelabel = " ".join( [ str(x) for x in self.diff_set(self.bags[x],self.bags[y] )] )
            out.write( "\tbag{} -> bag{}  [label=\"{}\"]\n".format(x+1,y+1,edgelabel) )

        out.write("\n}\n")

    ## @brief expand tree decomposition to have a certain maximal diff
    ## size
    ##
    ## the diff set is the set of variables in the child that do not
    ## occur in the parent bag
    ##
    ## @pre the tree is connected
    def expand_treedecomposition(self, maxdiffsize=1):
    
        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in range(0, len(l), n):
                yield l[i:i + n]
        
        n = len(self.bags)
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
                diff = self.diff_set(self.bags[u],self.bags[v])
                if len(diff) > maxdiffsize:
                    sep = self.sep_set(self.bags[u],self.bags[v])

                    if (u,v) in self.edges:
                        self.edges.remove((u,v))

                    if (v,u) in self.edges:
                        self.edges.remove((v,u))

                    last_bag_idx = u

                    newbag = sep
                    
                    for ext in chunks(diff[:-1],maxdiffsize):
                        newbag.extend(ext)
                        self.bags.append(newbag[:])
                        self.edges.append([ last_bag_idx, next_bag_idx])
                        last_bag_idx = next_bag_idx
                        next_bag_idx += 1

                    self.edges.append([ last_bag_idx, v])

        self.update_adjacency_lists()

############################################################
## writing graphs to files in different formats
##

def dotfile_to_pdf(graphfile):
    """
    Plot graph (written to a pdf file) by calling graphviz's dot tool
    on a given dot file
    """
    outfile = re.sub(r".dot$",".pdf",graphfile)
    subprocess.check_output(["dot","-Tpdf","-o",outfile,graphfile])

def write_dgf(out, num_nodes, edges):
    """
    Write graph in dgf format
    @param out output file handle
    @param num_nodes number of nodes
    @param edges the edges of the graph
    """

    edge_num=len(edges)

    out.write("p tw {} {}\n".format(num_nodes, edge_num))
    for (u,v) in sorted(edges):
        out.write("e {} {}\n".format(u+1,v+1))

def write_dot(out, num_nodes, edges):
    """
    Write graph in dot format
    @param out output file handle
    @param num_nodes number of nodes
    @param edges the edges of the graph
    """

    edge_num=len(edges)

    out.write("graph G{\n\n")

    for v in range(num_nodes):
        out.write("\tnode{idx} [label=\"{idx}\"]\n".format(idx=v))

    for (u,v) in edges:
        out.write("\tnode{} -- node{}\n".format(u,v))

    out.write("\n}\n")


############################################################
## Interface TDlib
##
## specific functions to call TDlib's tree decomposition tool
## and parse the result
##

def parseTD_TDlib(tdfh, num_nodes):
    """
    Parse tree decomposition as written by TDlib

    @param tdfh file handle to tree decomposition in dot format
    @param num_nodes number of nodes
    @returns tree decomposition

    assume td format as written by the TDlib program
    """

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

    return TreeDecomp(bags,edges)

## make tree decomposition from file using TDlib
def makeTDFile_TDlib( filename, num_nodes, edges,
                      *, strategy=2 ):
    """
    Compute tree decomposition of a graph by TDlib

    Generates tree decomposition and writes the result to file.

    @param filename base name for input .dfg and output .td files
    @param num_nodes number of nodes
    @param edges specifies edges of a graph; nodes are indexed 0-based
    @param strategy integer code of decomposition strategy:
      1) Permutation To Tree Decomposition with GreedyDegree
      2) Permutation To Tree Decomposition with GreedyFillIn
      3) Permutation To Tree Decomposition with Lex-BFS: Lexicographic Breadth First Search
      4) Permutation To Tree Decomposition with MCS; Maximum Cardinality Search Algorithm
      5) PermutationGuesser
    """

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

def makeTD_TDlib(num_nodes, edges, *, strategy=2):
    """
    Make tree decomposition using TDlib

    @todo generate unique / thread-safe tmp name
    """
    tmpfile = "/tmp/treedecomp-tmp~"
    makeTDFile_TDlib( tmpfile, num_nodes, edges, strategy=strategy )
    tdfile = tmpfile+".td"
    with open(tdfile) as tdfh:
        td =  parseTD_TDlib( tdfh, num_nodes )
    #os.remove(tdfile)
    return td

## End of TDlib-specific functions
##----------------------------------------------------------

############################################################
## Interface the htd library
##
## specific functions to interface libhtd
##

def makeTD_htd(num_nodes, edges, maxdiffsize=1):
    """
    Obtain tree decomposition using htd
    @param num_nodes number of nodes
    @param edges specifies edges of a graph; nodes are indexed 0-based
    """

    myhtd = htd.HTD(num_nodes,edges)
    myhtd.decompose()
    td = TreeDecomp(myhtd.bags(), myhtd.edges())

    td.expand_treedecomposition(maxdiffsize)

    return td

## End of libhtd-specific functions
##----------------------------------------------------------

############################################################
## Interface tree demposition libraries

def makeTD(num_nodes, edges, *, method=0, **kwargs):

    """
    Obtain tree decomposition

    Dispatches to tree decomp libraries based on 'method'

    @param num_nodes number of nodes
    @param edges specifies edges of a graph; nodes are indexed 0-based
    @param method chooses tree decomposition method (0=libhtd; othw. use TDlib with strategy=method
    @return bags,edges
    """

    if str(method) == "0":
        return makeTD_htd(num_nodes, edges, **kwargs)
    else:
        return makeTD_TDlib(num_nodes, edges, strategy=method)

## End of Interface tree demposition libraries
##----------------------------------------------------------

if __name__ == "__main__":
    pass
