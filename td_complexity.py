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


# Location of TDlib
tdlib_base = "."

import re
import argparse
import timeit
import os

import treedecomp as td 
import rnastuff as rna

str_to_dep = { "basepair": rna.structure_to_basepair_dependencies,
               "stacking": rna.structure_to_stacking_dependencies
}

def process_instance(infh,keep_graphs,plot_graphs,strategy,models):
    structures = rna.read_inp(infh)

    tw=dict()
    time=dict()

    for model in models:
        edges = list()
        for s in structures:
            sarray = rna.parseRNAStructure(s)
            str_to_dep[model](sarray,edges)
            
        seqlen = len(structures[0])

        edges = rna.unique_edges(edges)

        filename=model

        start_time = timeit.default_timer()
        tdfilename = td.makeTDFile(edges, filename, strategy=strategy, keep_graphs=keep_graphs, plot_graphs=plot_graphs)
        time[model] = timeit.default_timer() - start_time

        with open(tdfilename,"r") as tdfh:
            bags,edges = td.parseTD(seqlen,tdfh)
            tw[model] = td.treewidth(bags)
            with open(tdfilename+".dot","w") as out:
                td.writeTD(out,bags,edges)

        if plot_graphs:
            td.dotfile_to_pdf(tdfilename+".dot")
            
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
