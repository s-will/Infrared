#!/usr/bin/env python3

# -----------------------------
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

###############################
## @file
## Romy based on InfraRed
##
## Samples over alignments

import random
import argparse
import itertools
import os

import infrared as ir
import treedecomp
import rna_support as rna

import RNA

## @brief InfraRed function to control hamming distance
class HammingDistance(ir.Function):
    def __init__(self, i, j, weight):
        super().__init__([i, j])
        self.i = i
        self.j = j
        self.weight = weight
    def __call__(self, a): #overrides
        a = a.values()
        if (a[self.i]==a[self.j]):
            return self.weight
        else:
            return 1

## @brief GC content feature
class GCFeature(ir.Feature):
    def __init__(self, weight, target, tolerance):
        super().__init__( "GC", weight, target, tolerance)
    def eval(self, sample):
        return rna.GC_content("".join(sample)) * 100

## @brief Turner energy feature
class EnergyFeature(ir.Feature):
    def __init__(self, index, structure, weight, target, tolerance):
        super().__init__( ("E",index), weight, target, tolerance )
        self.structure = structure
        self.index = index

    def eval(self, sample):
        import RNA
        return RNA.energy_of_struct(sample[self.index], self.structure)

    def idstring(self):
        return "".join(map(str,self.identifier))

## @brief sequence distance feature
class DistanceFeature(ir.Feature):
    def __init__(self, phylotree_edge, weight, target, tolerance):
        i = phylotree_edge[0]
        j = phylotree_edge[1]

        super().__init__( ("D", i, j), weight, target, tolerance )

        self.i = i
        self.j = j

    @staticmethod
    def hamming_distance(xs,ys):
        return sum(map(lambda x:x[0]!=x[1], zip(xs,ys)))

    def eval(self, sample):
        return self.hamming_distance( sample[self.i], sample[self.j] )

    def idstring(self):
        return "_".join(map(str,self.identifier))


class RomyConstraintNetwork(ir.ConstraintNetwork):
    def __init__(self, seqnum, seqlen, phylotree, structure, features):

        super().__init__()

        self.seqnum = seqnum
        self.seqlen = seqlen
        self.phylotree = phylotree
        self.structure_string = structure
        self.structure = rna.parseRNAStructureBps(structure)
        self.features = features

        self.generate_cn()

    ## @brief Get variable id
    ## @param seq_id sequence id
    ## @param pos position
    def vid(self, seq_id, pos):
        return seq_id * self.seqlen + pos

    ## @brief functions for gc-content control
    def gc_control_functions(self):
        return 
    def generate_cn(self):
        # determine constraints due to consensus structure
        self.bp_constraints = [ ( [ self.vid(i,p), self.vid(i,q) ],
                                  [ ir.ComplConstraint(self.vid(i,p), self.vid(i,q)) ] )
                                for p, q in self.structure for i in range(self.seqnum) ]

        # set up dependencies due to evolutionary distance
        self.distance_functions = [ ( [ self.vid(i,k), self.vid(j,k)  ],
                                      [HammingDistance(self.vid(i,k), self.vid(j,k),
                                                       self.features[("D", i, j)].weight )] )
                                    for k in range(self.seqlen) for (i,j) in self.phylotree ] 
        
        # determine energy functions due to consensus structure
        self.energy_functions = [ ( [ self.vid(i,p), self.vid(i,q) ],
                                    [ir.BPEnergy(self.vid(i,p), self.vid(i,q),
                                                 not (p-1,q+1) in self.structure,
                                                 self.features[("E",i)].weight )] )
                                  for p,q in self.structure for i in range(self.seqnum) ]

        # GC content control
        self.gc_functions = [ ( [ self.vid(i,j) ],
                                [ ir.GCControl( self.vid(i,j), self.features["GC"].weight ) ] )
                              for i in range(self.seqnum) for j in range(self.seqlen) ]

        self.dependencies = list(map(lambda x: x[0],
                                     self.distance_functions
                                     + self.energy_functions))

class RomyTreeDecomposition(ir.TreeDecomposition):
    def __init__(self, cn, *, method):
        super().__init__(cn.seqlen*cn.seqnum, cn.dependencies, method=method)

        self.domains = 4
        self.cn = cn

    def get_bag_assignments(self):
        bagconstraints = self.assign_to_all_bags(self.cn.bp_constraints)

        functions = (self.cn.gc_functions
                     + self.cn.energy_functions
                     + self.cn.distance_functions)

        bagfunctions = self.assign_to_bags(functions)

        return bagconstraints, bagfunctions

class RomySampler(ir.MultiDimensionalBoltzmannSampler):
    def __init__( self, seqnum, seqlen, phylotree, structure, features, *, method ):
        super().__init__( features )

        self.seqnum = seqnum
        self.seqlen = seqlen
        self.phylotree = phylotree
        self.structure = structure
        self.features = features

        self.method = method

        self.setup_engine()

    ## @brief Generate constraint network
    ## @param features dictionary of features
    ## @return constraint network
    def gen_constraint_network(self, features):
        return RomyConstraintNetwork(self.seqnum, self.seqlen, self.phylotree,
                                     self.structure, self.features)

    ## @brief Generate tree decomposition
    ## @param cn constraint network
    ## @return tree decomposition
    def gen_tree_decomposition(self, cn):
        ## make tree decomposition
        return RomyTreeDecomposition( cn, method = self.method)

    def values_to_alignment(self,vals):
        seqs = []
        for i in range(self.seqnum):
            seq = ir.values_to_sequence(vals[:self.seqlen])
            vals = vals[self.seqlen:]
            seqs.append(seq)
        return seqs

    ## @brief Calculate sample
    ## @return sampled RNA sequence
    def sample(self):
        return self.values_to_alignment(super().sample().values())


## @brief command line tool definition
## @param args command line arguments
def main(args):

    ## init seed
    if args.seed == None:
        ir.seed(random.randint(0,2**31))
    else:
        ir.seed(args.seed)

    # set base pair energies
    ir.set_bpenergy_table( ir.params_bp )

    ## hard code an instance
    structure = "((((.((((...))))))))"
    seqnum = 5
    seqlen = len(structure)
    phylotree = [(3,0),(3,4),(4,1),(4,2)]

    features = [ GCFeature( 1, 66, 5 ) ] 
    features.extend( [ EnergyFeature( i, structure, 1, -10, 5 )
                       for i in range(seqnum) ] )
    features.extend( [ DistanceFeature( edge, 1, 2, 1 )
                       for edge in phylotree ] )    
    # transform features to the corresponding feature dictionary
    features = { f.identifier:f for f in features }

    sampler = RomySampler( seqnum, seqlen, phylotree, structure, features, method=args.method )

    ## optionally, write tree decomposition
    if args.plot_td:
        sampler.plot_td("treedecomp.dot")

    if args.verbose:
        print("Treewidth:",sampler.treewidth())


    fstats = ir.FeatureStatistics()

    ## sample

    if args.mdbs:
        sample_generator = sampler.targeted_samples()
    else:
        sample_generator = sampler.samples()

    sample_count = 0
    for alignment in sample_generator:
        print(alignment,end='')
        for i in range(seqnum):
            feat_id, value = fstats.record( sampler.features[("E",i)], alignment )
            print(" {}={:3.2f}".format(feat_id, value), end='')
 
        for (i,j) in sorted(phylotree):
            feat_id, value = fstats.record( sampler.features[("D",i,j)], alignment )
            print(" {}={:3.2f}".format(feat_id, value), end='')
            
        feat_id, value = fstats.record( sampler.features["GC"], alignment )
        print(" GC={:3.2f}".format(value),end='')

        print()

        sample_count += 1
        if sample_count >= args.number:
            break
 
    if args.verbose:
        if not fstats.empty():
            print("----------")
            print("Summary: ",end='')
        fstats.report()


if __name__ == "__main__":
    ## command line argument parser
    parser = argparse.ArgumentParser(description='Boltzmann sampling of homologous sequences')

    # parser.add_argument('infile', help="Input file")

    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; otherwise pass to TDlib as strategy)")

    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")

    parser.add_argument('--seed', type=int, default=None,
                        help="Seed infrared's random number generator (def=auto)")
    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")

    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")

    parser.add_argument('--mdbs', action="store_true",
                        help="Perform multi-dim Boltzmann sampling to aim at targets")

    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")

    args=parser.parse_args()

    main(args)
