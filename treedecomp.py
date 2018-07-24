#!/usr/bin/env python3

#
# Module to handle tree decompositions
# generate tree decompositions by calling TDlib
# plot dependency graphs and tree decompositions using dot
#
# @todo replace external program calls by using Python modules

tdlib_base = "."


import os
import subprocess
import re
from math import sqrt,ceil

def dotfile_to_pdf(graphfile):
    outfile = re.sub(r".dot$",".pdf",graphfile)
    subprocess.check_output(["dot","-Tpdf","-o",outfile,graphfile])

def parseTD(n, td):
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
    
    for i in range(n):
        if not i in present:
            bags.append([i])

    return (bags,edges)

def writeTD(out,bags,edges,offset=0):
    """
    Write tree decomposition in dot format
    @param out output file handle
    @param bags the bags of the TD
    @param edges the edges of the DT
    """

    def baglabel(bag):
        lwidth = ceil( sqrt(len(bag)) )
        lnum   = ceil( len(bag) / lwidth )
        xs = [str(i+offset) for i in bag]
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

def bagdifference(xs,ys):
    """Compute directed difference between bags ys-xs"""
    xdict = { x:1 for x in xs }
    return [ y for y in ys if not y in xdict ]

def treewidth(bags):
    """Determine tree width
    @param bags the bags of the tree decomposition
    """
    return max([len(bag) for bag in bags]) - 1

def write_dgf(edges,out):
    edge_num=len(edges)

    node_num=0
    if edge_num>0:
        node_num=max([ max(u,v)+1 for (u,v) in edges ])

    out.write("p tw {} {}\n".format(node_num, edge_num))
    for (u,v) in sorted(edges):
        out.write("e {} {}\n".format(u+1,v+1))

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

def inc_edges(edges):
    """Make 0-based edges 1-based"""
    return [ (i+1,j+1) for (i,j) in edges ]

def makeTDFile( edges, filename, *, strategy=2, keep_graphs=False, plot_graphs=False ):
    """Compute tree decomposition of a graph by TDlib
    
    Generates tree decomposition and writes the result to file.
    
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
               str(strategy), inname, outname
               #, filename+".dot"
        ]
        subprocess.check_output(cmd)

        if keep_graphs or plot_graphs:
            with open(dotfilename, "w") as dot:
                write_dot(inc_edges(edges), dot)

        if plot_graphs:
            dotfile_to_pdf(dotfilename)

        if plot_graphs and not keep_graphs:
            os.remove(dotfilename)
                
        os.remove(inname)

        return outname

def makeTD(size, edges, *, strategy=2):
    tmpfile = "/tmp/treedecomp-tmp~"
    makeTDFile(edges,tmpfile, strategy=strategy)
    tdfile = tmpfile+".td"
    with open(tdfile) as td:
        td =  parseTD(size,td)
    #os.remove(tdfile) 
    return td
        
if __name__ == "__main__":
    pass
