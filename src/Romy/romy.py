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

import numpy as np

import clustering as cl
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
from Bio import AlignIO, Phylo
import sys

import RNA

def hamming_distance(xs,ys):
    return sum(map(lambda x:x[0]!=x[1], zip(xs,ys)))

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
        fc = RNA.fold_compound(sample[self.index])
        return fc.eval_structure(self.structure)

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

    def eval(self, sample):
        return hamming_distance( sample[self.i], sample[self.j] )

    def idstring(self):
        return "_".join(map(str,self.identifier))

## @brief Construct and hold constraint network for mult-target design
## based in the base pair energy model
class RomyConstraintNetworkFactory:
    ## @brief Constructor
    def __init__( self ):
        pass

    ## @brief create constraint network
    ## @param seqlen length of sequences
    ## @param structures list of target structures in dot bracket format
    ## @param features the features containing weights of energies and GC control
    def create( self, seqsize, seqnum, seqlen, phylotree, structure, features ):

        self.seqsize = seqsize
        self.seqnum = seqnum
        self.seqlen = seqlen
        self.phylotree = phylotree
        self.structure_string = structure
        self.structure = rna.parseRNAStructureBps(structure)
        self.features = features

        self.generate_constraints_and_functions()

        cn = ir.ConstraintNetwork( varnum = seqlen * seqnum, domains = 4,
                constraints = self.bp_constraints,
                functions = self.gc_functions + self.energy_functions + self.distance_functions )

        return cn


    ## @brief Get variable id
    ## @param seq_id sequence id
    ## @param pos position
    def vid(self, seq_id, pos):
        return seq_id * self.seqlen + pos

    ## @brief Generate constraint network constraints and functions
    def generate_constraints_and_functions(self):
        # determine constraints due to consensus structure
        self.bp_constraints = [ rna.ComplConstraint(self.vid(i,p), self.vid(i,q))
                                for p, q in self.structure for i in range(self.seqnum) ]

        # set up dependencies due to evolutionary distance
        self.distance_functions = [ HammingDistance(self.vid(i,k), self.vid(j,k),
                                                    self.features[("D", i, j)].weight )
                                    for k in range(self.seqlen) for (i,j) in self.phylotree ]

        # determine energy functions due to consensus structure
        self.energy_functions = [ rna.BPEnergy(self.vid(i,p), self.vid(i,q),
                                               not (p-1,q+1) in self.structure,
                                               self.features[("E",i)].weight )
                                  for p,q in self.structure for i in range(self.seqsize) ]

        # GC content control
        self.gc_functions = [ rna.GCControl( self.vid(i,j), self.features["GC"].weight )
                              for i in range(self.seqnum) for j in range(self.seqlen) ]


class RomySampler(ir.MultiDimensionalBoltzmannSampler):
    def __init__( self, seqsize, seqnum, seqlen, phylotree, structure, features,
                  *, td_factory, cn_factory ):
        super().__init__( features )

        self.seqsize = seqsize
        self.seqnum = seqnum
        self.seqlen = seqlen
        self.phylotree = phylotree
        self.structure = structure
        self.features = features

        self.td_factory = td_factory
        self.cn_factory = cn_factory

        self.setup_engine()

    ## @brief Generate constraint network
    ## @param features dictionary of features
    ## @return constraint network
    def gen_constraint_network(self, features):
        return self.cn_factory.create(self.seqsize, self.seqnum, self.seqlen, self.phylotree,
                                     self.structure, self.features)

    ## @brief Generate tree decomposition
    ## @param cn constraint network
    ## @return tree decomposition
    def gen_tree_decomposition(self, cn):
        ## make tree decomposition
        return td_factory.create( cn )

    def values_to_alignment(self,vals):
        seqs = []
        for i in range(self.seqnum):
            seq = rna.values_to_sequence(vals[:self.seqlen])
            vals = vals[self.seqlen:]
            seqs.append(seq)
        return seqs

    def gen_cluster_tree( self ):
        return ir.ClusterTree( self.cn, td = self.td )

    ## @brief Calculate sample
    ## @return sampled RNA sequence
    def sample(self):
        return self.values_to_alignment(super().sample().values())

##Additional useful functions
def all_edges(tree,n):
    phylo_v = []
    phylo = []
    seqnum = n
    for clade in tree.find_clades(order='level'):
        for child in clade:
            if child.name[:5]=="Inner":
                if child.name=="Inner": #Case of only one Inner node
                    child_name=n
                    seqnum +=1
                else:
                    i = int(child.name.split('r')[1])+n-1
                    if i+1>seqnum:
                        seqnum=i+1
                    child_name = i
            else: child_name=int(child.name)

            if clade.name[:5]=="Inner":
                if clade.name=="Inner": #Case of only one Inner node
                    clade_name=n
                    seqnum +=1
                else:
                    i = int(clade.name.split('r')[1])+n-1
                    if i+1>seqnum:
                        seqnum=i+1
                    clade_name = i
            else: clade_name=int(clade.name)

            phylo_v.append((clade_name,child_name,child.branch_length))
            phylo.append((clade_name,child_name))


    return phylo_v,phylo,seqnum

def all_edges_nhx(tree,n):
    phylo_v = []
    phylo = []
    index_in, index_out = n,0
    for clade in tree.find_clades(order='level'):
        for child in clade:

            if clade.comment == None:
                clade_name = index_in
                index_in += 1
            else:
                clade_name = index_out
                index_out +=1

            if child.comment == None:
                child_name = index_in
                index_in += 1
            else:
                child_name = index_out
                index_out +=1

            phylo_v.append((clade_name,child_name,child.branch_length))
            phylo.append((clade_name,child_name))

    return phylo_v,phylo,index_in

def get_alignment_features(args):

    #aln=AlignIO.read(args.infile,'stockholm')
    sequences= list(RNA.file_msa_read(args.infile)[2])
    msa_size= len(sequences)
    if args.newick!=None:
        tree=Phylo.read(args.newick,'newick')
        phylo_v,phylotree,seqnum = all_edges_nhx(tree,msa_size)
    else:
        ds_mat = [[0 for i in range(i+1)] for i in range(msa_size)]
        for i in range(msa_size):
            for j in range(i):
                ds_mat[i][j] = hamming_distance(sequences[i],sequences[j])

        distance_matrix = DistanceMatrix([str(i) for i in range(msa_size)],ds_mat)

        constructor = DistanceTreeConstructor()
        tree=constructor.nj(distance_matrix)
        phylo_v,phylotree,seqnum = all_edges(tree,msa_size)





    #GC content and energy

    if args.struct == None:
        target_struct = RNA.alifold(sequences)[0]
    else:
        target_struct = args.struct
    #target_struct = "(((((((.((((.......))))((((((.......))))))...(((((.......))))))))))))."
    gc,energies=cl.analyze_alignments(sequences,target_struct)
    return {"Sequences": sequences, "Structure":target_struct,"GC":gc,"Energies":energies,"Tree":tree,"Phylotree":phylotree,"Phylo_v":phylo_v,"Seqnum":seqnum,"Size":msa_size}

def sequence_to_gap_pattern( ali_sequence ):
    """
    @brief convert a sequence to a gap pattern
    @param ali_sequence alignment sequence, possibly containing gaps '-'
    """
    return [ x=='-' for x in ali_sequence ]

"""
A quite specialized tree class

Supports traversal based on a list of edges
@note makes quite a few assumptions on edges and nodes (@see infer_inner_gap_patterns)!
"""
class Tree():
    def __init__(self,edges):
        self.nodes = set()
        for (i,j) in edges:
            self.nodes.add(i)
            self.nodes.add(j)

        self.adjacency = dict()
        for (i,j) in edges:
            if i not in self.adjacency:
                 self.adjacency[i]=list()
            if j not in self.adjacency:
                self.adjacency[j]=list()
            self.adjacency[i].append(j)
            self.adjacency[j].append(i)

        self.parent = dict()
        self.children = dict()

        self.root = max(self.nodes)
        self._init_parents_and_children(self.root)

    def _init_parents_and_children(self,i,parent=None):
        self.parent[i] = parent
        self.children[i] = []
        for j in self.adjacency[i]:
            if j!=parent:
                self._init_parents_and_children(j,i)
                self.children[i].append(j)

def infer_inner_gap_patterns( sequences, tree_edges ):
    """
    @brief Infer the gap patterns at the inner leaves of the tree
    @param sequences sequences/alignment strings of the alignment (including gaps)
    @param tree as list of edges
    @returns list of all gap patterns

    @pre the indices of sequences and node indices of the leave nodes in tree correspond;
    in tree the leaves must have the indices in range(number of leaves), inner nodes have integer indices
    in range(number of leaves, number of nodes); all sequences have the same length
    """
    print("infer_inner_gap_patterns",tree_edges)

    leave_gap_patterns = [ sequence_to_gap_pattern(x) for x in sequences ]
    leave_num = len(leave_gap_patterns)

    if leave_num==0:
        return []

    seqlen = len(sequences[0])

    tree = Tree(tree_edges)

    ## We run a fitch maximum parsimony algo with trace back, separately on each alignment column
    ##

    # we follow this schema:

    # foreach alignment column
    #   ## fitch_fwd
    #   traverse nodes in post order, at each node
    #      determine and store max parsimonius scores for each node type
    #
    #   determine best node type
    #
    #   ## fitch_tb
    #   traverse pre-order, at each node
    #      choose types for children that yield maximum score

    values = [False,True]

    cost_tab = dict()

    def fitch_fwd( node, col ):
        ncost_tab = dict()
        if tree.children[node] == []:
            # init the table
            g = leave_gap_patterns[node][col]
            ncost_tab[g] = 0
            ncost_tab[not g] = 10**9 # hackish for 'infinite cost'
        else:
            # run algo on kids and infer the table
            for child in tree.children[node]:
                fitch_fwd(child, col)

            for v in values:
                ncost_tab[v] = sum( min( cost_tab[child][v], cost_tab[child][not v] + 1 ) for child in tree.children[node] )

        cost_tab[node] = ncost_tab

    tb_values = dict()
    def fitch_tb( node, val ):
        tb_values[node]=val

        # pick optimal values for children
        children = tree.children[node]

        if children == []:
            return

        best_cost = 10**9
        best_values = None
        for children_values in itertools.product( values, repeat=len(children) ):
            cost = sum( cost_tab[children[cidx]][children_values[cidx]] + (children_values[cidx]!=val) for cidx in range(len(children)) )
            if cost < best_cost:
                best_cost = cost
                best_values = children_values

        for cidx in range(len(children)):
            fitch_tb( children[cidx], best_values[cidx] )

    gap_patterns = [ [] for node in tree.nodes ]

    for col in range(seqlen):
        fitch_fwd( tree.root, col )

        best_cost,best_val = min( (cost_tab[tree.root][val],val) for val in values )

        fitch_tb( tree.root, best_val )

        for node in tree.nodes:
            gap_patterns[node].append(tb_values[node])

    return gap_patterns

def gap_pattern_to_string(gap_pattern):
    def f(x):
        return {True:'-',False:'.'}[x]
    return "".join( [ f(x) for x in gap_pattern ] )

def subtract_gap_distances( tree, gps ):
    """
    @brief subtract hamming distances between gap patterns for each branch of a tree

    Corrects the tree lengths to reflect hamming distances between non-gap positions

    @param tree tree edges with lengths (list of triples x,y,l)
    @param gps gap patterns of the nodes
    """
    return [ (x,y,l - hamming_distance(gps[x],gps[y])) for (x,y,l) in tree ]

## @brief command line tool definition
## @param args command line arguments
def main(args):

    if args.listtds:
        print("Avalaible tree decomposition methods", treedecomp.get_td_factory_descriptors())
        return

    ## init seed
    if args.seed == None:
        ir.seed(random.randint(0,2**31))
    else:
        ir.seed(args.seed)

    # set base pair energies
    rna.set_bpenergy_table()

    ## hard code an instance

    """     structure = "((((.((((...))))))))"
    seqnum = 5
    seqlen = len(structure)
    phylotree = [(3,0),(3,4),(4,1),(4,2)] #List of edges

    features = [ GCFeature( 1, 66, 5 ) ] #GC content of 66% with a tolerance of 5%
    features.extend( [ EnergyFeature( i, structure, 1, -10, 5 ) #Energy of -10 with a tolerance of 5%
                       for i in range(seqnum) ] )
    features.extend( [ DistanceFeature( edge, 1, 2, 1 ) #We want to have a hamming distance of 2% for each sequence
                       for edge in phylotree ] )  """

    #Get instance from alignment
    msa_features = get_alignment_features(args)

    sequences = msa_features["Sequences"]
    structure=msa_features["Structure"]

    seqnum=msa_features["Seqnum"]
    seqsize=msa_features["Size"]
    seqlen=len(structure)
    phylotree=msa_features["Phylotree"]
    phylotree_v=msa_features["Phylo_v"]

    # infer gap patterns at inner nodes
    gap_patterns = infer_inner_gap_patterns(sequences, phylotree)
    if args.verbose:
        print([ gap_pattern_to_string(x) for x in gap_patterns ])

    # update the branch lengths in the phylo tree by subtracting gap-caused distances
    corrected_phylotree_v = subtract_gap_distances( phylotree_v, gap_patterns )
    if args.verbose:
        print(phylotree_v,'--->',corrected_phylotree_v)

    ## GC feature
    features = [ GCFeature(args.gc_weight,msa_features["GC"],args.gc_tolerance) ]

    ## Energy features
    features.extend( [ EnergyFeature( i, structure, args.energy_weight, msa_features["Energies"][i], args.energy_tolerance )
                       for i in range(msa_features["Size"]) ] )

    ## Distance features
    features.extend( [ DistanceFeature((edge[0],edge[1]), args.distance_weight, edge[2], args.distance_tolerance )
                       for edge in corrected_phylotree_v ] )


    # transform features to the corresponding feature dictionary
    features = { f.identifier:f for f in features }

    # handle args.td
    td_factory = treedecomp.td_factory_from_descriptor(args.td)
    if td_factory is None:
        sys.stderr.write("[ERROR] Invalid tree decomposition method: "+args.td+"\n")
        exit(-1)

    cn_factory = RomyConstraintNetworkFactory()

    sampler = RomySampler( seqsize, seqnum, seqlen, phylotree, structure, features,
                           td_factory = td_factory,
                           cn_factory = cn_factory
                           )

    ## optionally, write tree decomposition
    if args.plot_td:
        sampler.plot_td("treedecomp.dot")

    if args.verbose:
        print("Treewidth:",sampler.treewidth())
        print("The targeted structure is: ",structure)
        print("The targeted average GC content is: ",msa_features["GC"])
        print("The targeted energies are: ",msa_features["Energies"])
    """         print("The phylogenetic tree will be printed in another window, close it tp have the results")
        Phylo.draw(msa_features["Tree"]) """


    fstats = ir.FeatureStatistics()

    ## sample

    if args.mdbs:
        sample_generator = sampler.targeted_samples()
    else:
        sample_generator = sampler.samples()

    sample_count = 0
    alignments=[]
    for alignment in sample_generator:
        msa = alignment[:seqsize]
        alignments.append(msa)
        print(msa,end='')
        for i in range(seqsize):
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

    return alignments,seqnum

if __name__ == "__main__":
    ## command line argument parser
    parser = argparse.ArgumentParser(description='Boltzmann sampling of homologous sequences')

    parser.add_argument('infile',type=str, help="Input Stockholm file of the alignment")

    parser.add_argument('--struct', type=str, default=None, help="Consensus structure for the alignment")

    parser.add_argument('--gc_tolerance', type=float, default=5, help="Target tolerance for the GC content")

    parser.add_argument('--energy_tolerance', type=float, default=5, help="Target tolerance for energies")

    parser.add_argument('--distance_tolerance', type=float, default=1, help="Target tolerance for hamming distances")

    parser.add_argument('--gc_weight', type=float, default=1, help="GC weight")

    parser.add_argument('--energy_weight', type=float, default=1, help="Energy weight")

    parser.add_argument('--distance_weight', type=float, default=1, help="Distance weight")

    parser.add_argument('--td', type=str, default="nx",
                        help="Method for tree decomposition (see --listtds)")
    parser.add_argument('--listtds', action="store_true", help="List available tree decomposition methods")

    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")

    parser.add_argument('--seed', type=int, default=None,
                        help="Seed infrared's random number generator (def=auto)")


    parser.add_argument("--newick",type=str,default=None,
                        help="Filename of the newick phylogenetic tree to use")

    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")

    parser.add_argument('--mdbs', action="store_true",
                        help="Perform multi-dim Boltzmann sampling to aim at targets")

    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")

    args=parser.parse_args()

    main(args)
