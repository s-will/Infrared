#!/usr/bin/env python3

# -----------------------------
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#
# Redprint provides Boltzmann sampling of sequences targeting multiple RNA structures.
#

###############################
## @file
## Redprint v2 based on InfraRed
##
## Defines the redprint library module and the redprint command line
## tool
##
## @note Dependencies: the redprint tool needs ViennaRNA's RNA pyhton
## module, which could require to set the python path like: export
## PYTHONPATH=$VRNAVIENNA_HOME/lib/python3.6/site-packages. For
## further dependencies see Infrared/treedecomp.py.

import random
import argparse
import itertools
import os

import libinfrared as ir
import treedecomp
import rna_support as rna


## Parameters for the base pair model (magic params from the Redprint paper)
params_bp={
    "GC_IN": -2.10208,
    "AU_IN": -0.52309,
    "GU_IN": -0.88474,
    "GC_TERM": -0.09070,
    "AU_TERM": 1.26630,
    "GU_TERM": 0.78566
}

## Parameters for the stacking model (magic params from the Redprint paper)
params_stacking = {
    "AUAU": -0.18826,
    "AUCG": -1.13291,
    "AUGC": -1.09787,
    "AUGU": -0.38606,
    "AUUA": -0.26510,
    "AUUG": -0.62086,
    "CGAU": -1.11752,
    "CGCG": -2.23740,
    "CGGC": -1.89434,
    "CGGU": -1.22942,
    "CGUA": -1.10548,
    "CGUG": -1.44085,
    "GUAU": -0.55066,
    "GUCG": -1.26209,
    "GUGC": -1.58478,
    "GUGU": -0.72185,
    "GUUA": -0.49625,
    "GUUG": -0.68876
}

############################################################
## Redprint Library

## @brief Convert integer (variable value) to nucleotide
## @note encoding A=0, C=1, G=2, U=3
def val2nucl(x):
    return "ACGU"[x]

## @brief Convert list of integers (variable values) to string
## (sequence) of nucleotides
def values2seq(xs):
    return "".join(map(val2nucl, xs))

## @brief set the bp energy table for Infrared
def set_bpenergy_table(tab):
    ir.BPEnergy.set_energy_table(tab)

## @brief set the stacking energy table for Infrared
def set_stacking_energy_table(tab):
    ir.StackEnergy.set_energy_table(tab)

## @brief A constraint network for multi-target RNA design
##
## 'Abstract' class, which provides some common functionality for
## specific energy models
##
## The constraint network specifies the number of variables, the
## specific dependencies, the constraints and functions of a
## multi-target RNA design instance
class RNAConstraintNetwork:
    def __init__(self, seqlen, structures, weights, gcweight):
        self.seqlen = seqlen
        self.structures = list(structures)
        self.weights = list(weights)
        self.gcweight = gcweight

        # to be filled later
        self.bpdependencies = []
        self.constraints=[]
        self.functions=[]

    ## @brief Generate base pair dependencies
    ##
    ## @return unique binary dependencies due to base pairs as list of lists of indices;
    def gen_bp_dependencies(self):
        bps = set()
        for structure in self.structures:
            for bp in structure:
                bps.add(bp)
        return [ [i,j] for (i,j) in list(bps) ]

    ## @brief Generate stacking dependencies
    ##
    ## @return unique 4-ary dependencies due to stacks of base pairs as list of lists of indices
    def stacking_dependencies(self):
        stacked_bps = set()
        for structure in self.structures:
            structureset = set(structure)
            for (i,j) in structure:
                if (i+1,j-1) in structureset:
                    stacked_bps.add((i,j))
        return [ [i,j,i+1,j-1] for (i,j) in list(stacked_bps) ]


    ## @brief Compute complementarity classes
    ##
    ## Compute complementarity classes for each connected component,
    ## i.e. the components of the bipartitiion induced by the
    ## complementarity constraints on base pairs
    ##
    ## Returns result in self.compl_classes
    def compute_compl_classes(self):
        ## Transform pairs into dict with keys=first and value=list of second components
        def accumulate_dict(xys):
            d = {x:[] for x,y in xys}
            for (x,y) in xys:
                d[x].append(y)
            return d

        ## perform depth-first traversal
        self.compl_classes = dict()

        visited = set()

        other_ends = accumulate_dict(self.bpdependencies
                                     + list( map(list,map(reversed,self.bpdependencies)) ))
        color = 1

        for x in range(self.seqlen):
            if x in visited: continue

            ## new component, color it
            self.compl_classes[x] = color
            color+=1

            if x not in other_ends:
                continue # if there is no base pair

            stack = [x]

            while stack:
                x = stack.pop()
                if x in visited: continue
                visited.add(x)
                for y in other_ends[x]:
                    self.compl_classes[y]= -self.compl_classes[x]
                    stack.append(y)

    ## @brief Add the functions for gc-content control to self.functions
    def add_gc_control(self):
        gc_control_funs = [ ( [i], [ir.GCControl(i,self.gcweight)] ) for i in range(self.seqlen) ]
        self.functions.extend( gc_control_funs )

    ## @brief Add the complementarity constraints to self.constraints
    def compl_constraints(self):
        # generate constraints and functions; assign them to bags
        return [ ([i,j], [ir.ComplConstraint(i,j)]) for [i,j] in self.bpdependencies ]

    ## @brief Removes redundant dependencies
    ##
    ## @param deps list of dependencies (where a dependency is a list of indices)
    ##
    ## Removes all dependencies that are subsumed by other dependencies in deps
    ##
    ## @return pruned list of dependencies
    @staticmethod
    def remove_redundant_dependencies(deps):
        def sublist(xs,ys):
            return all(x in ys for x in xs)
        def subsumed(dep,deps):
            return any( len(dep)<len(dep2) and sublist(dep,dep2) for dep2 in deps )
        return [ dep for dep in deps if not subsumed(dep,deps) ]

## @brief Construct and hold constraint network for mult-target design based in the base pair energy model
class RNAConstraintNetworkBasePair(RNAConstraintNetwork):
    ## @brief Constructor from sequence length, structures and weights
    ##
    ## @param seqlen length of sequences
    ## @param structures list of target structures in dot bracket format
    ## @param weights list of weights for the single structures
    ## @param gcweight weight for GC-content
    def __init__(self, seqlen, structures, weights, gcweight):

        super().__init__(seqlen, structures, weights, gcweight)

        self.generate_cn_basepair_model()

        self.compute_compl_classes()

    ## @brief Generate constraint network for the base pair model
    def generate_cn_basepair_model(self):
        self.bpdependencies = self.gen_bp_dependencies()
        self.dependencies = self.bpdependencies

        self.constraints = self.compl_constraints()

        self.functions = list()

        for k,structure in enumerate(self.structures):
            structureset = set(structure)
            self.functions.extend( [ ( [i,j],
                                       [ir.BPEnergy(i,j, not (i-1,j+1) in structureset, self.weights[k])] )
                                     for (i,j) in structure ] )

        self.add_gc_control()

## @brief Construct and hold constraint network for mult-target design based in the stacking energy model
class RNAConstraintNetworkStacking(RNAConstraintNetwork):
    ## @brief Constructor from sequence length, structures and weights
    ##
    ## @param seqlen length of sequences
    ## @param structures list of target structures in dot bracket format
    ## @param weights list of weights for the single structures
    ## @param gcweight weight for GC-content
    def __init__(self, seqlen, structures, weights, gcweight):

        super().__init__(seqlen, structures, weights, gcweight)

        self.generate_cn_stacking_model()

        self.compute_compl_classes()

    ## @brief Generate constraint network for the stacking model
    def generate_cn_stacking_model(self):
        self.bpdependencies = self.gen_bp_dependencies()
        self.dependencies = self.remove_redundant_dependencies( self.stacking_dependencies() + self.bpdependencies )


        self.constraints = self.compl_constraints()

        self.functions = list()

        for k,structure in enumerate(self.structures):
            structureset = set(structure)
            self.functions.extend( [ ( [i,j,i+1,j-1],
                                       [ir.StackEnergy(i, j, self.weights[k])] )
                                     for (i,j) in structure
                                     if (i+1,j-1) in structureset
            ] )

        self.add_gc_control()

## @brief Holds tree decomposition for RNA design, constructs RNA design cluster tree
##
class RNATreeDecomposition:

    ## @brief Constructor from constraint network
    ##
    ## @param cn the RNA design constraint network
    ## @param add_redundant_constraints whether to insert redundant
    ## constraints (in practice, one usually wants this for performance!)
    ## @param method tree decomposition method
    ##
    ## Calls tree decomposition algorithm (according to method)
    def __init__(self, cn, *, add_redundant_constraints=True, method=0):
        self.cn = cn

        self.add_redundant_constraints = add_redundant_constraints

        # from dependencies generate list of binary edges
        bindependencies  = self.expand_to_cliques(cn.dependencies)

        # generate tree decomposition -> bags, edges  (beware: translate between 0/1-based)
        self.td = treedecomp.makeTD(cn.seqlen, bindependencies, method = method)

        self.bagsets = list(map(set,self.td.bags))

    ## @brief Expand non-binary dependencies to cliques of binary deps
    ## @param dependencies list of dependencies
    ## @return list of binary dependencies
    @staticmethod
    def expand_to_cliques(dependencies):
        bindeps = list()
        for d in dependencies:
            bindeps.extend( itertools.combinations(d,2) )
        return bindeps


    ## @brief Get the indices of all bags that contain a set of variables
    ## @param bvars the set of variables
    ## @return list of indices of the bags that contain bvars
    def find_all_bags(self,bvars):
        return [ i for i,bag in enumerate(self.bagsets) if all( x in bag for x in bvars ) ]

    ## @brief Find a bag that contains a set of variables
    ## @param bvars the set of variables
    ## @return index of first bag that contains bvars (or None if there is none)
    def find_bag(self,bvars):
        bags = self.find_all_bags(bvars)
        if len(bags)>0:
            return bags[0]
        else:
            return None

    ## @brief assign constraints or functions to bags
    ## @param constraints list of constraints
    ## @return list where constraints are placed at corresponding bag indices
    ##
    ## Assigns such that each constraint/function is assigned to exactly one bag
    ## that contains its dependencies
    ##
    ## @pre for each constraint/function there is one bag that contains its dependencies;
    ## otherwise the constraint/function is not assigned
    def assign_to_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bagsets)) }
        for (cvars,ccons) in constraints:
            bagconstraints[self.find_bag(cvars)].extend(ccons)
        return bagconstraints

    ## @brief assign constraints or functions to all possible bags
    ## @param constraints list of constraints
    ## @return list where constraints are placed at corresponding bag indices
    ##
    ## Assigns constraint/function is to each bags that contains its
    ## dependencies.
    def assign_to_all_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bagsets)) }
        for (cvars,ccons) in constraints:
            for bidx in self.find_all_bags(cvars):
                bagconstraints[bidx].extend(ccons)
        return bagconstraints

    ## @brief Fillin complementarity class constraints
    ##
    ## Fills in constraints due to complementarity to all bags. Adds
    ## SameComplClass or DifferentComplClass constraints, wherever a
    ## variable would have to be enumerated unnconstrained otherwise
    ## @param classes complementarity classes @param
    ## existing_constraints @param[in,out] bagconstraints
    ##
    ## Avoids to add constraints which already exist
    def fillin_class_constraints(self, classes, existing_constraints, bagconstraints):
        existing = set( (vars[0],vars[1]) for (vars,constr) in existing_constraints )
        for bagidx,bag in enumerate(self.bagsets):
            bagvars=sorted(list(bag))
            for j in bagvars[1:]:
                if all( (i,j) not in existing for i in range(0,j) ):
                    i = bagvars[0]
                    if classes[i] == classes[j]:
                        #print("Add same:",i,j)
                        bagconstraints[bagidx].append(ir.SameComplClassConstraint(i,j))
                    elif classes[i] == - classes[j]:
                        #print("Add diff:",i,j)
                        bagconstraints[bagidx].append(ir.DifferentComplClassConstraint(i,j))

    ## @brief Construct the cluster tree
    def construct_cluster_tree(self):
        if self.add_redundant_constraints:
            bagconstraints = self.assign_to_all_bags(self.cn.constraints)
            self.fillin_class_constraints(self.cn.compl_classes,
                                          self.cn.constraints,
                                          bagconstraints)
        else:
            bagconstraints = self.assign_to_bags(self.cn.constraints)

        bagfunctions = self.assign_to_bags(self.cn.functions)

        ct = ir.ClusterTree(self.cn.seqlen, 4);

        children = set() # keep record of all non-root nodes (which have been seen as children)
        for bagidx in self.td.toposorted_bag_indices():
            if not bagidx in children:
                # enumerate subtree
                stack = [(None,bagidx)]
                while stack:
                    (p,i) = stack[-1]
                    stack = stack[:-1]
                    bagvars = sorted(list(self.bagsets[i]))

                    if p==None:
                        cluster = ct.add_root_cluster(bagvars)
                    else:
                        cluster = ct.add_child_cluster(p,bagvars)

                    for x in bagconstraints[i]: #  + self.cn.gen_redundant_bpconstraints(bagvars)
                        ct.add_constraint(cluster, x)

                    for x in bagfunctions[i]:
                        ct.add_function(cluster, x)

                    for j in self.td.adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct


# END Redprint Library
# ##########################################################

## @brief Redprint command line tool definition
## @param args command line arguments
def main(args):
    import RNA

    ## init seed
    if args.seed == None:
        ir.seed(random.randint(0,2**31))
    else:
        ir.seed(args.seed)

    # set base pair energies ( AU,GC,GU in stems and terminal )
    set_bpenergy_table(
        list(map(lambda x: params_bp[x],
                 [ "AU_IN", "GC_IN", "GU_IN",
                   "AU_TERM", "GC_TERM", "GU_TERM" ] )))

    set_stacking_energy_table(
        list(map(lambda x: params_stacking[x],
                 [ "AUAU", "AUUA",
                   "AUCG", "AUGC",
                   "AUGU", "AUUG",

                   "CGAU", "CGUA",
                   "CGCG", "CGGC",
                   "CGGU", "CGUG",

                   "GUAU", "GUUA",
                   "GUCG", "GUGC",
                   "GUGU", "GUUG" ] )))

    ## read instance
    with open(args.infile) as infh:
        structures = rna.read_inp(infh)

    weights=args.weight
    if weights is None: weights=[1]

    if len(weights) < len(structures):
        weights.extend([weights[-1]]*(len(structures)-len(weights)))

    if len(structures) == 0:
        print("At least one structure is required in input")
        exit(-1)

    seqlen = len(structures[0])

    #parse structures
    bps = list(map(rna.parseRNAStructureBps,structures))

    ## build constraint network
    if args.model in ["bp","basepair"]:
        cn = RNAConstraintNetworkBasePair( seqlen, bps,
                                           weights, args.gcweight )
    elif args.model in ["stack","stacking"]:
        cn = RNAConstraintNetworkStacking( seqlen, bps,
                                           weights, args.gcweight )
    else:
        print("Model",args.model,"unknown! Please see help for supported modules.")
        exit(-1)

    #print(sorted([ vars for (vars,c) in cn.constraints]))

    ## make tree decomposition
    rtd = RNATreeDecomposition( cn,
                               add_redundant_constraints = not args.no_redundant_constraints,
                               method=args.method )

    #print(rtd.td.bags)

    ## optionally, write tree decomposition
    if args.plot_td:
        dotfilename = "treedecomp.dot"
        with open(dotfilename,"w") as dot:
            rtd.td.writeTD(dot)
        treedecomp.dotfile_to_pdf(dotfilename)
        os.remove(dotfilename)

    if args.verbose:
        print("Treewidth:",rtd.td.treewidth())

    ## make cluster tree
    ct = rtd.construct_cluster_tree()

    ## evaluate
    ct.evaluate()

    # for statistics
    features = dict()
    def register_feature(feat_id,val):
        if feat_id not in features:
            features[feat_id] = []
        features[feat_id].append(val)

    ## sample
    for x in range(0,args.number):
        sample = ct.sample()
        seq = values2seq(sample.values())
        print(seq,end='')
        if args.turner:
            for i,struc in enumerate(structures):
                eos = RNA.energy_of_struct(seq,struc)
                feat_id = "E_{}".format(i+1)
                register_feature(feat_id,eos)
                print(" {}={:3.2f}".format(feat_id,eos),end='')

        if args.gc:
            gc = rna.GC_content(seq)
            feat_id = "GC"
            register_feature(feat_id,gc)
            print(" GC={:3.2f}".format(gc*100),end='')

        if args.checkvalid:
            for i,struc in enumerate(structures):
                if not rna.is_valid(seq, struc):
                    print(" INVALID{}: {}".format(i,str(rna.invalid_bps(seq,struc))),end='')
        print()

    if args.verbose:
        if features:
            print("----------")
            print("Summary: ",end='')

        def mean(xs):
            xs = list(xs)
            return sum(xs)/len(xs)

        for fid in features:
            mu = mean(features[fid])
            std = (mean(map(lambda x: x**2,features[fid])) - mu**2 )**0.5

            print(" {}={:3.2f} +/-{:3.2f}".format(fid,mu,std),end='')
        print()

if __name__ == "__main__":
    ## Redprint command line argument parser
    parser = argparse.ArgumentParser(description='Boltzmann sampling for RNA design with multiple target structures')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; otherwise pass to TDlib as strategy)")
    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")
    parser.add_argument('--seed', type=int, default=None, help="Seed infrared's random number generator (def=auto)")

    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    parser.add_argument('--model', type=str, default="bp",
                        help="Energy model used for sampling [bp=base pair model,stack=stacking model]")
    parser.add_argument('--turner', action="store_true",
                        help="Report Turner energies of the single structures for each sample")
    parser.add_argument('--gc', action="store_true",
                        help="Report GC contents of the single structures for each sample")

    parser.add_argument('--checkvalid', action="store_true",
                        help="Check base pair complementarity for each structure and for each sample")

    parser.add_argument('--no_redundant_constraints', action="store_true", help="Do not add redundant constraints")

    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")
    parser.add_argument('-w','--weight', type=float, action="append", help="Structure weight (def=1; in case, use last weight for remaining structures)")

    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")


    args=parser.parse_args()

    main(args)
