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


import random
import argparse
from collections import Counter
import itertools
import os

import libinfrared as ir
import treedecomp
import rna_support as rna
import treedecomp

# importing Vienna package if not in path could require setting python path like
# export PYTHONPATH=$HOME/Soft/ViennaRNA-2.4.8/lib/python3.6/site-packages
import RNA

##############################
# The magic parameters used in the redprint paper
#
GC_IN = -2.10208
AU_IN = -0.52309
GU_IN = -0.88474
GC_TERM = -0.09070
AU_TERM = 1.26630
GU_TERM = 0.78566
#
STACK_AUAU = -0.18826
STACK_AUCG = -1.13291
STACK_AUGC = -1.09787
STACK_AUGU = -0.38606
STACK_AUUA = -0.26510
STACK_AUUG = -0.62086
STACK_CGAU = -1.11752
STACK_CGCG = -2.23740
STACK_CGGC = -1.89434
STACK_CGGU = -1.22942
STACK_CGUA = -1.10548
STACK_CGUG = -1.44085
STACK_GUAU = -0.55066
STACK_GUCG = -1.26209
STACK_GUGC = -1.58478
STACK_GUGU = -0.72185
STACK_GUUA = -0.49625
STACK_GUUG = -0.68876


############################################################
## Redprint Library

def val2nucl(x):
    return "ACGU"[x]
def values2seq(xs):
    return "".join(map(val2nucl,xs))

def set_bpenergy_table(tab):
    ir.BPEnergy.set_energy_table(tab)

def set_stacking_energy_table(tab):
    ir.StackEnergy.set_energy_table(tab)

def accumulate_dict(xys):
    d = { x:[] for x,y in xys }
    for (x,y) in xys:
        d[x].append(y)
    return d

class RNAConstraintNetwork:
    def __init__(self, seqlen, structures, weights, gcweight):
        self.seqlen = seqlen
        self.structures = list(structures)
        self.weights = list(weights)
        self.gcweight = gcweight

    ## generate the dependencies due to base pairs
    ## (makes dependencies unique)
    def bp_dependencies(self):
        bps = set()
        for structure in self.structures:
            for bp in structure:
                bps.add(bp)
        return [ [i,j] for (i,j) in list(bps) ]

    ## generate the dependencies due to base pairs
    ## (makes dependencies unique)
    def stacking_dependencies(self):
        stacked_bps = set()
        for structure in self.structures:
            structureset = set(structure)
            for (i,j) in structure:
                if (i+1,j-1) in structureset:
                    stacked_bps.add((i,j))
        return [ [i,j,i+1,j-1] for (i,j) in list(stacked_bps) ]


    ## compute complementarity classes for each connected component,
    ## i.e. the components of the bipartitiion induced by the
    ## complementarity constraints on base pairs
    ##
    ## @writes result in self.compl_classes
    def compute_compl_classes(self):
        ## perform depth-first traversal
        self.compl_classes = dict()

        visited = set()

        other_ends = accumulate_dict(self.bpdependencies + list( map(list,map(reversed,self.bpdependencies)) ))

        color = 1
        
        for x in range(self.seqlen):
            if x in visited: continue
            if x not in other_ends: continue # if there is no base pair

            ## new component, color it
            stack = [x]
            self.compl_classes[x] = color
            color+=1

            while stack:
                x = stack.pop()
                if x in visited: continue
                visited.add(x)
                for y in other_ends[x]:
                    self.compl_classes[y]= -self.compl_classes[x]
                    stack.append(y)
    
    # add the functions for gc-content control to self.functions 
    def add_gc_control(self):
        gc_control_funs = [ ( [i], [ir.GCControl(i,self.gcweight)] ) for i in range(self.seqlen) ]
        self.functions.extend( gc_control_funs )
        
    # add the complementarity constraints to self.constraints
    def compl_constraints(self):
        # generate constraints and functions; assign them to bags
        return [ ([i,j], [ir.ComplConstraint(i,j)]) for [i,j] in self.bpdependencies ]

    @staticmethod
    def remove_redundant_dependencies(deps):
        def sublist(xs,ys):
            return all(x in ys for x in xs)
        def subsumed(dep,deps):
            return any( len(dep)<len(dep2) and sublist(dep,dep2) for dep2 in deps )
        return [ dep for dep in deps if not subsumed(dep,deps) ]


class RNAConstraintNetworkBasePair(RNAConstraintNetwork):
    def __init__(self, seqlen, structures, weights, gcweight):

        super().__init__(seqlen, structures, weights, gcweight)

        self.generate_cn_basepair_model()
        
        self.compute_compl_classes()
        

    ## generate constraint network for the base pair model
    def generate_cn_basepair_model(self):
        self.bpdependencies = self.bp_dependencies()
        self.dependencies = self.bpdependencies

        self.constraints = self.compl_constraints()

        self.functions = list()
        
        for k,structure in enumerate(self.structures):
            structureset = set(structure)
            self.functions.extend( [ ( [i,j],
                                       [ir.BPEnergy(i,j, not (i-1,j+1) in structureset, self.weights[k])] ) 
                                     for (i,j) in structure ] )

        self.add_gc_control()


class RNAConstraintNetworkStacking(RNAConstraintNetwork):
    def __init__(self, seqlen, structures, weights, gcweight):

        super().__init__(seqlen, structures, weights, gcweight)

        self.generate_cn_stacking_model()
        
        self.compute_compl_classes()
        

    ## generate constraint network for the base pair model
    def generate_cn_stacking_model(self):
        self.bpdependencies = self.bp_dependencies()
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

class RNATreeDecomposition:
    def __init__(self, cn, *, add_redundant_constraints=True, method):
        self.cn = cn

        self.add_redundant_constraints = add_redundant_constraints

        # from dependencies generate list of binary edges
        bindependencies  = self.expand_to_cliques(cn.dependencies)

        # generate tree decomposition -> bags, edges  (beware: translate between 0/1-based)
        bags,edges = treedecomp.makeTD(cn.seqlen, bindependencies, method = method)

        self.edges = edges
        self.bags = list(map(set,bags))

    @staticmethod
    def expand_to_cliques(dependencies):
        bindeps = list()
        for d in dependencies:
            bindeps.extend( itertools.combinations(d,2) )
        return bindeps

    
    ## get (first) index of bag that contains all variables
    def find_all_bags(self,bvars):
        return [ i for i,bag in enumerate(self.bags) if all( x in bag for x in bvars ) ]

    def find_bag(self,bvars):
        bags = self.find_all_bags(bvars)
        if len(bags)>0:
            return bags[0]
        else:
            return None

    # asisgn constraints or functions to bags
    # @returns list where constraints are placed at corresponding bag indices
    def assign_to_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bags)) }
        for (cvars,ccons) in constraints:
            bagconstraints[self.find_bag(cvars)].extend(ccons)
        return bagconstraints

    def assign_to_all_bags(self,constraints):
        bagconstraints = { i:[]  for i in range(len(self.bags)) }
        for (cvars,ccons) in constraints:
            for bidx in self.find_all_bags(cvars):
                bagconstraints[bidx].extend(ccons)
        return bagconstraints

    ## Fills in constraints due to complementarity to all bags. Adds
    ## SameComplClass or DifferentComplClass constraints, wherever a
    ## variable would have to be enumerated unnconstrained otherwise
    ## @param classes complementarity classes
    ## @param existing_constraints
    ## @param[in,out] bagconstraints
    ##
    ## avoids to add constraints which already exist
    def fillin_class_constraints(self, classes, existing_constraints, bagconstraints):
        existing = set( (vars[0],vars[1]) for (vars,constr) in existing_constraints )
        for bagidx,bag in enumerate(self.bags):
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

    @staticmethod
    def toposort(n,adj):
        visited = set()
        sorted = list()

        def toposort_helper(i):
            visited.add(i)
            for j in adj[i]:
                if not j in visited:
                    toposort_helper(j)
            sorted.append(i)

        for i in range(n):
            if not i in visited:
                toposort_helper(i)
        return sorted[::-1]

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

        ## perform topological sort
        adj = { i:[] for i in range(len(self.bags))}
        for [i,j] in self.edges:
            adj[i].append(j)

        sorted_bags = self.toposort(len(self.bags),adj)

        children = set() # keep record of all non-root nodes (which have been seen as children)
        for bagidx in sorted_bags:
            if not bagidx in children:
                # enumerate subtree
                stack = [(None,bagidx)]
                while stack:
                    (p,i) = stack[-1]
                    stack = stack[:-1]
                    bagvars = sorted(list(self.bags[i]))

                    if p==None:
                        cluster = ct.add_root_cluster(bagvars)
                    else:
                        cluster = ct.add_child_cluster(p,bagvars)

                    for xcon in bagconstraints[i]: #  + self.cn.gen_redundant_bpconstraints(bagvars)
                        ct.add_constraint(cluster, xcon)

                    for xfun in bagfunctions[i]:
                        ct.add_function(cluster, xfun)

                    #print(self.bags[i],list(map(lambda x:x.vars(),bagconstraints[i])))

                    for j in adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct


## END Redprint Library
############################################################

def main(args):
    ## init seed
    ir.seed(random.randint(0,2**31))

    # set base pair energies ( AU,GC,GU in stems and terminal )
    set_bpenergy_table( [ AU_IN,
                              GC_IN,
                              GU_IN,
                              AU_TERM,
                              GC_TERM,
                              GU_TERM ] )

    set_stacking_energy_table( [ STACK_AUAU, STACK_AUUA,
                                     STACK_AUCG, STACK_AUGC,
                                     STACK_AUGU, STACK_AUUG,
                                     
                                     STACK_CGAU, STACK_CGUA,
                                     STACK_CGCG, STACK_CGGC,
                                     STACK_CGGU, STACK_CGUG,
                                     
                                     STACK_GUAU, STACK_GUUA,
                                     STACK_GUCG, STACK_GUGC,
                                     STACK_GUGU, STACK_GUUG ] )

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
    td = RNATreeDecomposition( cn, 
                               add_redundant_constraints = not args.no_redundant_constraints,
                               method=args.method )

    #print(td.bags)

    ## optionally, write tree decomposition
    if args.plot_td:
        dotfilename = "treedecomp.dot"
        with open(dotfilename,"w") as dot:
            treedecomp.writeTD(dot, td.bags, td.edges)
        treedecomp.dotfile_to_pdf(dotfilename)
        os.remove(dotfilename)


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
        seq = values2seq(sample.values())
        print(seq,end='')
        if args.turner:
            for i,struc in enumerate(structures):
                print(" E_{}={:3.2f}".format(i+1,RNA.energy_of_struct(seq,struc)),end='')

        if args.checkvalid:
            for i,struc in enumerate(structures):
                if not rna.is_valid(seq, struc):
                    print(" INVALID{}: {}".format(i,str(rna.invalid_bps(seq,struc))),end='')

        print()

        ## statistics
        for i in range(seqlen):
            counters[i][seq[i]] += 1

    if args.verbose:
        print(sum(counters[1:],counters[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate training data')
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; othwerwise pass to TDlib as strategy)")
    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")
    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    parser.add_argument('--model', type=str, default="bp",
                        help="Energy model used for sampling [bp=base pair model,stack=stacking model]")
    parser.add_argument('--turner', action="store_true",
                        help="Report Turner energies of the single structures for each sample")
    parser.add_argument('--checkvalid', action="store_true",
                        help="Check base pair complementarity for each structure and for each sample")

    parser.add_argument('--no_redundant_constraints', action="store_true", help="Do not add redundant constraints")

    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")
    parser.add_argument('-w','--weight', type=float, action="append", help="Structure weight (def=1; in case, use last weight for remaining structures)")

    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")


    args=parser.parse_args()

    main(args)
