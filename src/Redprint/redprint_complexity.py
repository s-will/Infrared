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


import re
import argparse
import timeit
import os

import random
from libinfrared import seed

import treedecomp as td 
import rna_support as rna

str_to_dep = { "basepair": rna.structure_to_basepair_dependencies,
               "stacking": rna.structure_to_stacking_dependencies
}

def process_instance(infh,keep_graphs,plot_graphs,method,models):
    structures = rna.read_inp(infh)

    tw=dict()
    time=dict()

    for model in models:
        edges = list()
        for s in structures:
            sarray = rna.parseRNAStructure(s)
            str_to_dep[model](sarray,edges)
            
        seqlen = len(structures[0])

        in_edges = rna.unique_edges(edges)

        filename=model

        ## optionally, write dependeny graph
        if keep_graphs or plot_graphs:
            dotfilename = filename+".dot"
            with open(dotfilename, "w") as dot:
                td.write_dot(dot, seqlen, in_edges)
            if plot_graphs:
                td.dotfile_to_pdf(dotfilename)
            if not keep_graphs:
                os.remove(dotfilename)
        
        ## compute decomposition and measure run time
        start_time = timeit.default_timer()
        bags,edges = td.makeTD(seqlen, in_edges,
                               method=method)
        time[model] = timeit.default_timer() - start_time

        ## record tree width
        tw[model] = td.treewidth(bags)

        ## optionally, write tree decomposition
        if keep_graphs or plot_graphs:
            dotfilename = filename+"_td.dot"
            with open(dotfilename,"w") as dot:
                td.writeTD(dot, bags, edges)
            if plot_graphs:
                td.dotfile_to_pdf(dotfilename)
            if not keep_graphs:
                os.remove(dotfilename)
        
    return  [ str(tw[model]) for model in models ] + [ str(time[model]) for model in models ] ;

if __name__ == "__main__":
    ## init seed
    seed(random.randint(0,2**31))

    parser = argparse.ArgumentParser(description='Generate training data')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--keep_graphs', action="store_true",
                        help="Keep the generated graphs and tree decompositions")
    parser.add_argument('--plot_graphs', action="store_true",
                        help="Plot generated dependency graphs and tree decompositions")
    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; othwerwise pass to TDlib as strategy)")

    args=parser.parse_args()

    models = str_to_dep.keys()

    with open(args.infile) as infh:
        results = process_instance(infh,args.keep_graphs,args.plot_graphs,args.method,models)
        print("\t".join([args.infile]+results))
