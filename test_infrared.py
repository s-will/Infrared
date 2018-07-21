#!/usr/bin/env python3

import infrared_module as irm

import argparse
from collections import Counter
import rna

def main(args):
    ## init seed
    irm.seed(1)

    INF = 1e6
    energy_tab = [INF,INF,INF, -2,  #A
                  INF,INF, -3,INF,  #C
                  INF, -3,INF, -1,  #G
                   -2,INF, -1,INF] #U
    
    irm.set_bpenergy_table(energy_tab)

    ## read instance
    with open(args.infile) as infh:
        structures = rna.read_inp(infh)

    if len(args.weight) != len(structures):
        print("Wrong number of weights specified for this instance, need",len(structures))
        exit(-1)

    if len(structures) == 0:
        print("At least one structure is required in input")
        exit(-1)

    seqlen = len(structures[0])

    # read model
    structures = map(rna.parseRNAStructureBps,structures)

    ## build constraint network
    cn = irm.RNAConstraintNetworkBasePair( seqlen, structures, args.weight, args.gcweight )

    ## make tree decomposition
    td = irm.RNATreeDecomposition( cn, strategy=args.strategy )

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
        seq = irm.values2seq(sample.values())
        print(seq)

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

    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")
    parser.add_argument('-w','--weight', type=float, action="append", help="Structure weight")

    args=parser.parse_args()

    main(args)
