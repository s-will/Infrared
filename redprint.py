#!/usr/bin/env python3

############################### 
# Redprint v2 based on InfraRed
##-----------------------------
# (C) Sebastian Will, 2018
# 
# This file is part of the InfraRed source code.
# 
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#
# Redprint provides Boltzmann sampling of sequences targeting multiple RNA structures.
# 
# This file defines the Redprint main function.
# 


import redprint_lib as rpl

import random
import argparse
from collections import Counter
from rna_support import read_inp, parseRNAStructureBps, is_valid, invalid_bps

# importing Vienna package if not in path could require setting python path like
# export PYTHONPATH=$HOME/Soft/ViennaRNA-2.4.8/lib/python3.6/site-packages
import RNA

def main(args):
    ## init seed
    rpl.seed(random.randint(0,2**31))

    energy_tab = [ -2, -3, -1 ]

    rpl.set_bpenergy_table(energy_tab)

    ## read instance
    with open(args.infile) as infh:
        structures = read_inp(infh)

    if args.weight is None or len(args.weight) != len(structures):
        print("Wrong number of weights specified for this instance, need",len(structures))
        exit(-1)

    if len(structures) == 0:
        print("At least one structure is required in input")
        exit(-1)

    seqlen = len(structures[0])

    #parse structures
    bps = list(map(parseRNAStructureBps,structures))

    ## build constraint network
    cn = rpl.RNAConstraintNetworkBasePair( seqlen, bps, 
                                           args.weight, args.gcweight )

    #print(sorted([ vars for (vars,c) in cn.constraints]))

    ## make tree decomposition
    td = rpl.RNATreeDecomposition( cn, 
                                   add_redundant_constraints = not args.no_redundant_constraints,
                                   strategy=args.strategy )

    #print(td.bags)

    treewidth = max(map(len,td.bags))-1
    print("Treewidth:",treewidth)
    
    ## make cluster tree
    ct = td.construct_cluster_tree()

    ## evaluate
    ct.evaluate()

    ## sample

    # for statistics
    counters=[Counter() for i in range(0,seqlen)]

    for x in range(0,args.number):
        sample = ct.sample()
        seq = rpl.values2seq(sample.values())
        print(seq,end='')
        if args.turner:
            for i,struc in enumerate(structures):
                print(" E_{}={:3.2f}".format(i+1,RNA.energy_of_struct(seq,struc)),end='')

        if args.checkvalid:
            for i,struc in enumerate(structures):
                if not is_valid(seq, struc):
                    print(" INVALID{}: {}".format(i,str(invalid_bps(seq,struc))),end='')

        print()

        ## statistics
        for i in range(seqlen):
            counters[i][seq[i]] += 1

    if args.verbose:
        print(sum(counters[1:],counters[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate training data')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--strategy', default="2", help="Strategy for tree decomposition (integer code; passed to TDlib)")
    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")
    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    parser.add_argument('--turner', action="store_true", help="Report Turner energies of the single structures for each sample")
    parser.add_argument('--checkvalid', action="store_true", help="Check base pair complementarity for each structure and for each sample")

    parser.add_argument('--no_redundant_constraints', action="store_true", help="Do not add redundant constraints")

    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")
    parser.add_argument('-w','--weight', type=float, action="append", help="Structure weight")

    args=parser.parse_args()

    main(args)
