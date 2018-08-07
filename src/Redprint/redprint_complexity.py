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

import libinfrared
import treedecomp
import rna_support as rna

str_to_dep = { "basepair": rna.structure_to_basepair_dependencies,
               "stacking": rna.structure_to_stacking_dependencies
}

def process_instance(infh,models,args):
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

        ## optionally, write dependeny graph
        if args.keep_graphs or args.plot_graphs:
            dotfilename = filename+".dot"
            with open(dotfilename, "w") as dot:
                treedecomp.write_dot(dot, seqlen, edges)
            if args.plot_graphs:
                treedecomp.dotfile_to_pdf(dotfilename)
            if not args.keep_graphs:
                os.remove(dotfilename)

        ## compute decomposition and measure run time
        start_time = timeit.default_timer()
        td = treedecomp.makeTD(seqlen, edges, method=args.method, maxdiffsize=args.maxdiffsize)
        time[model] = timeit.default_timer() - start_time

        ## record tree width0
        tw[model] = td.treewidth()

        ## optionally, write tree decomposition
        if args.keep_graphs or args.plot_graphs:
            dotfilename = filename+"_td.dot"
            with open(dotfilename,"w") as dot:
                td.writeTD(dot)
            if args.plot_graphs:
                treedecomp.dotfile_to_pdf(dotfilename)
            if not args.keep_graphs:
                os.remove(dotfilename)

    return  [ str(tw[model]) for model in models ] + [ str(time[model]) for model in models ] ;

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Report details on redprints complexity')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--keep_graphs', action="store_true",
                        help="Keep the generated graphs and tree decompositions")
    parser.add_argument('--plot_graphs', action="store_true",
                        help="Plot generated dependency graphs and tree decompositions")
    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; othwerwise pass to TDlib as strategy)")
    parser.add_argument('--seed', type=int, default=None, help="Seed infrared's random number generator. Concerns only method=0 (def=auto)")

    parser.add_argument('--maxdiffsize', type=int, default=1,
                        help="Maximum size of diff sets (Concerns only method=0)")

    parser.add_argument('--model', type=str, default="bp",
                        help="Energy model used for sampling [bp=base pair model,stack=stacking model]")

    args=parser.parse_args()

    ## init seed
    if args.seed == None:
        libinfrared.seed(random.randint(0,2**31))
    else:
        libinfrared.seed(args.seed)

    ## model
    if args.model is None:
        models = str_to_dep.keys()
    elif args.model in ["bp","basepair"]:
        models = ["basepair"]
    elif args.model in ["stack","stacking"]:
        models = ["stacking"]
    else:
        print("Model",args.model,"unknown! Please see help for supported modules.")
        exit(-1)

    with open(args.infile) as infh:
        results = process_instance(infh,models,args)
        print("\t".join([args.infile]+results))
