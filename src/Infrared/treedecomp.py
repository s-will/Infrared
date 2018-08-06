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
## writing graphs to files in different formats
##

def dotfile_to_pdf(graphfile):
    """
    Plot graph (written to a pdf file) by calling graphviz's dot tool
    on a given dot file
    """
    outfile = re.sub(r".dot$",".pdf",graphfile)
    subprocess.check_output(["dot","-Tpdf","-o",outfile,graphfile])

def writeTD(out, bags, edges, offset=0):
    """
    Write tree decomposition in dot format
    @param out output file handle
    @param bags the bags of the TD
    @param edges the edges of the DT
    """

    def baglabel(bag):
        if len(bag)==0:
            return ""
        
        lwidth = ceil( sqrt(len(bag)) )
        lnum   = ceil( len(bag) / lwidth )
        xs = [str(i+offset) for i in sorted(bag)]
        lines=list()
        for i in range(0,lnum):
            lines.append(" ".join(xs[i*lwidth:(i+1)*lwidth]))
        return "\\n".join(lines)

    out.write("graph G {\n\n")

    for bagid,bag in enumerate(bags):
        label = baglabel(bag)
        out.write( "\tbag{} [label=\"{}\"]\n".format(bagid+1,label) )

    out.write("\n\n")

    for (x,y) in edges:
        edgelabel = " ".join( [ str(x) for x in bagdifference(bags[x],bags[y] )] )
        out.write( "\tbag{} -- bag{}  [label=\"{}\"]\n".format(x+1,y+1,edgelabel) )

    out.write("\n}\n")

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

##############################
## convenience functions for the used representations of trees and
## tree decompositions
##

def max_node(edges):
    """
    Maximum node occuring in a list of edges
    """
    maxnode=0
    if len(edges)>0:
        maxnode=max([ max(u,v) for (u,v) in edges ])
    return maxnode

def bagdifference(xs,ys):
    """
    Compute directed difference between bags ys-xs
    """
    xdict = { x:1 for x in xs }
    return [ y for y in ys if not y in xdict ]

def treewidth(bags):
    """
    Determine tree width
    @param bags the bags of the tree decomposition
    """
    return max([len(bag) for bag in bags]) - 1

# def inc_edges(edges):
#     """Make 0-based edges 1-based"""
#     return [ (i+1,j+1) for (i,j) in edges ]
# 

############################################################
## Interface TDlib
##
## specific functions to call TDlib's tree decomposition tool
## and parse the result
##

def parseTD_TDlib(td, num_nodes):
    """
    Parse tree decomposition as written by TDlib

    @param td file handle to tree decomposition in dot format
    @param num_nodes number of nodes
    @returns tree decomposition, i.e. list of bags and edges between bags

    assume td format as written by the TDlib program
    """

    bags = list()
    edges = list()

    for line in td:
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

    return (bags,edges)

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
    with open(tdfile) as td:
        td =  parseTD_TDlib( td, num_nodes )
    #os.remove(tdfile)
    return td

## End of TDlib-specific functions
##----------------------------------------------------------

############################################################
## Interface the htd library
##
## specific functions to interface libhtd
##

def makeTD_htd(num_nodes, edges):
    """
    Obtain tree decomposition using htd
    @param num_nodes number of nodes
    @param edges specifies edges of a graph; nodes are indexed 0-based
    """

    myhtd = htd.HTD(num_nodes,edges)
    myhtd.decompose()

    return (myhtd.bags(), myhtd.edges())

## End of libhtd-specific functions
##----------------------------------------------------------

############################################################
## Interface tree demposition libraries

def makeTD(num_nodes, edges, *, method=0):

    """
    Obtain tree decomposition

    Dispatches to tree decomp libraries based on 'method'

    @param num_nodes number of nodes
    @param edges specifies edges of a graph; nodes are indexed 0-based
    @param method chooses tree decomposition method (0=libhtd; othw. use TDlib with strategy=method
    @return bags,edges
    """

    if str(method) == "0":
        return makeTD_htd(num_nodes, edges)
    else:
        return makeTD_TDlib(num_nodes, edges, strategy=method)

## End of Interface tree demposition libraries
##----------------------------------------------------------

if __name__ == "__main__":
    pass
