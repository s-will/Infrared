#!/usr/bin/env python3

#
# Determine the tree width for instances to multi-target structure design
#
# Author: Sebastian Will, Oct 2017
#
# Input is read from .inp files
# USAGE: td_complexity.py <inpfile>
#
# Determine tree width after decomp with TDlib in two energy models
#  -- (weighted) base pair energy model
#  -- stacking energy model
#
# The execution time of TDlib is logged.
#
# RESULT:
# writes tab-separated line
# filename tw_bp tw_st time_bp time_st
#
#
# INPUT FORMAT:
# * inp files contain the structures in dot bracket notation, one per line
# * structures must start with . or ( or )
# * pseudoknots are supported; possible bracket symbols: (),[],{},<>
# * file content after the first ; is ignored
#
# Side effects: (over)writes files basepair.* and stacking.* in current directory,
# which are removed on correct termination
# --> current directory MUST be writable
#
#
# Location of TDlib
tdlib_base = "../../RNARedPrint/lib"


import os
import subprocess
import re
import argparse
import timeit
from math import sqrt,ceil

def parseRNAStructure(structure):
    """Parse RNA structure including pseudoknots
    @param structure
    @returns array
    """

    opening = "([{<"
    closing = ")]}>"

    stack = { op:list() for op in opening }
    bps = [-1]*len(structure)

    for i,c in enumerate(structure):
        for (op,cl) in zip(opening,closing):
            if c==op:
                stack[op].append(i)
            elif c==cl:
                j = stack[op].pop()
                bps[i] = j
                bps[j] = i

    return bps

def unique_edges(xs):
    d = { (min(x,y),max(x,y)):1 for (x,y) in xs }
    return d.keys()

def structure_to_basepair_dependencies(structure,acc=[]):
    """Convert structure to list of edges in basepair model"""

    edges = acc

    for (i,j) in enumerate(structure):
        if i<j:
            edges.append((i,j))

    return edges

def structure_to_stacking_dependencies(structure,acc=[]):
    """Convert structure to list of edges in stacking model"""

    edges = acc

    for (i,j) in enumerate(structure):
        if i+1<j-1 and structure[i+1]==j-1:
            edges.extend([(i,j),(i+1,j-1),(i+1,j),(i,j-1),(i,i+1),(j,j-1)])

    return edges

str_to_dep = { "basepair": structure_to_basepair_dependencies,
               "stacking": structure_to_stacking_dependencies
}


def inc_edges(edges):
    """Make 0-based edges 1-based"""

    return [ (i+1,j+1) for (i,j) in edges ]

def write_dgf(edges,out):
    edge_num=len(edges)

    node_num=0
    if edge_num>0:
        node_num=max([ max(u,v) for (u,v) in edges ])

    out.write("p tw {} {}\n".format(node_num, edge_num))
    for (u,v) in sorted(edges):
        out.write("e {} {}\n".format(u,v))

def write_dot(edges,out):
    edge_num=len(edges)

    node_num=0
    if edge_num>0:
        node_num=max([ max(u,v) for (u,v) in edges ])

    out.write("graph G{\n\n")

    for v in range(node_num):
        out.write("\tnode{idx} [label=\"{idx}\"]\n".format(idx=v+1))
         
    for (u,v) in edges:
        out.write("\tnode{} -- node{}\n".format(u,v))

    out.write("\n}\n")

def plot_graph(graphfile):
    outfile = re.sub(r".dot$",".pdf",graphfile)
    subprocess.check_output(["dot","-Tpdf","-o",outfile,graphfile])

def bagdifference(xs,ys):
    """Compute directed difference between bags ys-xs"""
    xdict = { x:1 for x in xs }
    return [ y for y in ys if not y in xdict ]

def parseTD(td):
    """Determine tree width
    @param td file handle to tree decomposition in dot format
    @returns tree decomposition, i.e. list of bags and edges between bags

    assume td format as written by the TDlib program
    """

    bags = list()
    edges = list()

    for line in td:
        # ignore bags "New_vertex"
        if re.search("New_Vertex",line): continue

        m = re.search(r"bag(\d+).+label=\"(.+)\"",line)
        if m:
            bagid=int(m.group(1))
            label = m.group(2)
            labels = re.findall('\d+', label)
            labels = [ int(label)-1 for label in labels ]
            if bagid!=len(bags)+1:
                raise IOError("Bag indices in td file must be consecutive (at bag {})!".format(bagid))
            bags.append(labels)
        else:
            m = re.search(r"bag(\d+) -- bag(\d+)", line)
            if m:
                edges.append((int(m.group(1))-1,
                              int(m.group(2))-1))
    return (bags,edges)

def writeTD(out,bags,edges):
    """
    Write tree decomposition in dot format
    @param out output file handle
    @param bags the bags of the TD
    @param edges the edges of the DT
    """

    def baglabel(bag):
        lwidth = ceil( sqrt(len(bag)) )
        lnum   = ceil( len(bag) / lwidth )
        xs = [str(i+1) for i in bag]
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
        edgelabel = " ".join( [ str(x+1) for x in bagdifference(bags[x],bags[y] )] )
        out.write( "\tbag{} -- bag{}  [label=\"{}\"]\n".format(x+1,y+1,edgelabel) )

    out.write("\n}\n")

def treewidth(bags):
    """Determine tree width
    @param bags the bags of the tree decomposition
    """
    return max([len(bag) for bag in bags]) - 1

def makeTD(edges,filename,strategy,keep_graphs, plot_graphs):
    """Compute tree decomposition of a graph by TDlib

    @param edges specifies edges of a graph; nodes are indexed (1-based)
    @param filename base name for input .dfg and output .td files
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
        write_dgf(edges,dgf)
        dgf.close()

        cmd = ["java","-cp", tdlib_base+"/treewidth-java",
               "nl.uu.cs.treewidth.TreeDecomposer",
               strategy, inname, outname
               #, filename+".dot"
        ]
        #print(" ".join(cmd))
        subprocess.check_output(cmd)

        if keep_graphs or plot_graphs:
            with open(dotfilename, "w") as dot:
                write_dot(edges, dot)

        if plot_graphs:
            plot_graph(dotfilename)

        if plot_graphs and not keep_graphs:
            os.remove(dotfilename)
                
        os.remove(inname)

        return outname

def read_inp(inpfh):
    """read structures from inp file
    @param fh input file handle to inp file
    """

    structures=list()

    for line in inpfh:
        if re.match(re.compile("[\(\)\.]"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(";", line, flags=0):
            break

    return structures

def process_instance(infh,keep_graphs,plot_graphs,strategy,models):
    structures = read_inp(infh)

    tw=dict()
    time=dict()

    for model in models:
        edges = list()
        for s in structures:
            sarray = parseRNAStructure(s)
            str_to_dep[model](sarray,edges)

        edges = unique_edges( inc_edges(edges) )

        filename=model

        start_time = timeit.default_timer()
        tdfilename = makeTD(edges, filename, strategy, keep_graphs, plot_graphs)
        time[model] = timeit.default_timer() - start_time

        with open(tdfilename,"r") as td:
            bags,edges = parseTD(td)
            tw[model] = treewidth(bags)
            #print(bags,edges)
            with open(tdfilename+".dot","w") as out:
                writeTD(out,bags,edges)

        if plot_graphs:
            plot_graph(tdfilename+".dot")
            
        if not keep_graphs:
            os.remove(tdfilename)
            os.remove(tdfilename+".dot")

    return  [ str(tw[model]) for model in models ] + [ str(time[model]) for model in models ] ;

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate training data')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--keep_graphs', action="store_true",
                        help="Keep the generated graphs and tree decompositions")
    parser.add_argument('--plot_graphs', action="store_true",
                        help="Plot generated dependency graphs and tree decompositions")
    parser.add_argument('--strategy', default="2", help="Strategy for tree decomposition (integer code; passed to TDlib)")

    args=parser.parse_args()

    models = str_to_dep.keys()

    with open(args.infile) as infh:
        results = process_instance(infh,args.keep_graphs,args.plot_graphs,args.strategy,models)
        print("\t".join([args.infile]+results))
