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
## Redprint v2 based on InfraRed
##
## Redprint provides Boltzmann sampling of sequences targeting
## multiple RNA structures.  This file defines the redprint library
## module and the redprint command line tool.
##
## @note Dependencies: the redprint tool needs ViennaRNA's RNA pyhton
## module, which could require to set the python path like: export
## PYTHONPATH=$VRNAVIENNA_HOME/lib/python3.6/site-packages. For
## further dependencies see Infrared/treedecomp.py.

import random
import argparse
import itertools
import os

import infrared as ir
import rna_support as rna

############################################################
## Redprint Library


## @brief A constraint network for multi-target RNA design
##
## 'Abstract' class, which provides some common functionality for
## specific energy models
##
## The constraint network specifies the number of variables, the
## specific dependencies, the constraints and functions of a
## multi-target RNA design instance
class RNAConstraintNetwork(ir.ConstraintNetwork):
    ## @brief Constructor
    ## @param seqlen seduence length
    ## @param structures structures as lists of base pairs
    ## @param features the features containing weights of energies and GC control
    def __init__( self, seqlen, structures, features ):
        super().__init__()
	## sequence length
        self.seqlen = seqlen
	## list of structures as lists of base pairs
        self.structures = list(structures)
	## features containing energy weights and GC control
        self.features = features

        # to be filled later
	## list of base pair dependencies
        self.bpdependencies = []
	## list of #ired.Constraint
	## @todo document Constraint object
        self.constraints=[]
	## list of #ired.Function
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
    ## @return unique 4-ary dependencies due to stacks of base pairs
    ## as list of lists of indices
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
        # Transform pairs into dict with keys=first and value=list of second components
        def accumulate_dict(xys):
            d = {x:[] for x,y in xys}
            for (x,y) in xys:
                d[x].append(y)
            return d

	## @brief complementarity classes in dictionary
	##
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
        gc_control_funs = [ ( [i], [ ir.GCControl( i, self.features["GC"].weight ) ] )
                            for i in range(self.seqlen) ]
        self.functions.extend( gc_control_funs )

    ## @brief Add the complementarity constraint for each base pair dependenciy to self.constraints
    ## @see ired.rnadesign.ComplConstraint
    def compl_constraints(self):
        # generate constraints and functions; assign them to bags
        return [ ( [i,j], [ ir.ComplConstraint(i,j) ] ) for [i,j] in self.bpdependencies ]


## @brief Construct and hold constraint network for mult-target design
## based in the base pair energy model
class RNAConstraintNetworkBasePair(RNAConstraintNetwork):
    ## @brief Constructor
    ##
    ## @param seqlen length of sequences
    ## @param structures list of target structures in dot bracket format
    ## @param features the features containing weights of energies and GC control
    def __init__( self, seqlen, structures, features ):

        super().__init__( seqlen, structures, features )

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
            self.functions.extend( [ ( [ i, j ],
                                       [ ir.BPEnergy( i,j, not (i-1,j+1) in structureset,
                                                      self.features[("E",k)].weight ) ] )
                                     for (i,j) in structure ] )

        self.add_gc_control()

## @brief Construct and hold constraint network for mult-target design
## based in the stacking energy model
class RNAConstraintNetworkStacking(RNAConstraintNetwork):
    ## @brief Constructor
    ##
    ## @param seqlen length of sequences
    ## @param structures list of target structures in dot bracket format
    ## @param features the features containing weights of energies and GC control
    def __init__(self, seqlen, structures, features ):

        super().__init__( seqlen, structures, features )

        self.generate_cn_stacking_model()

        self.compute_compl_classes()

    ## @brief Generate constraint network for the stacking model
    def generate_cn_stacking_model( self ):
        self.bpdependencies = self.gen_bp_dependencies()
        self.dependencies = self.remove_redundant_dependencies( self.stacking_dependencies() + self.bpdependencies )


        self.constraints = self.compl_constraints()

        self.functions = list()

        for k,structure in enumerate(self.structures):
            structureset = set(structure)
            self.functions.extend( [ ( [ i,j,i+1,j-1 ],
                                       [ ir.StackEnergy( i, j, self.features[("E",k)].weight ) ] )
                                     for (i,j) in structure
                                     if (i+1,j-1) in structureset
            ] )

        self.add_gc_control()

## @brief Holds tree decomposition for RNA design, constructs RNA design cluster tree
##
class RNATreeDecomposition(ir.TreeDecomposition):

    ## @brief Constructor from constraint network
    ##
    ## @param cn the RNA design constraint network
    ## @param add_red_constrs whether to insert redundant
    ## constraints (in practice, one usually wants this for performance!)
    ## @param method tree decomposition method
    ##
    ## Calls tree decomposition algorithm (according to method)
    def __init__(self, cn, *, add_red_constrs=True, method=0):
        super().__init__(cn.seqlen, cn.dependencies, method=method)

        self.domains = 4

        self.cn = cn

        self.add_red_constrs = add_red_constrs

    ## @brief assign functions and constraints to bags
    ##
    ## @return bagconstraints, bagfuntions
    ##
    ## the function is called by construct_cluster_tree
    def get_bag_assignments(self):
        if self.add_red_constrs:
            bagconstraints = self.assign_to_all_bags(self.cn.constraints)
            self.fillin_class_constraints(self.cn.compl_classes,
                                          self.cn.constraints,
                                          bagconstraints)
        else:
            bagconstraints = self.assign_to_bags(self.cn.constraints)

        bagfunctions = self.assign_to_bags(self.cn.functions)

        return bagconstraints, bagfunctions

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


## @brief GC content feature
##
## Defines the feature 'GC content'
## @see infrared.Feature
class GCFeature(ir.Feature):

    ## The constructor.
    ## The identifier is set to "GC" by default
    def __init__(self, weight, target=None, tolerance=None):
        super().__init__( "GC", weight, target, tolerance)


    ## GC content of a given sample
    def eval(self, sample):
        return rna.GC_content(sample) * 100

## @brief Turner energy feature
##
## Defines the feature 'Turner energy'
##
## @note This feature exemplifies that the evaluation by the functions
## that control the sampling together with the weight of this feature
## does not have to be identical to the value of the feature: in
## RedPrint, the functions evaluate to some simplified energy in the
## base pair or stacking model, while the feature value is the Turner
## energy of the sequence as computed by RNA.energy_of_struct.
## @see infrared.Feature
class EnergyFeature(ir.Feature):
    ## The constructor
    ## @param index structure index
    ## @param structure secondary structure
    def __init__(self, index, structure, weight, target=None, tolerance=None):
        super().__init__( ("E",index), weight, target, tolerance )
        self.structure = structure
    ## Compute Turner energy of given sample
    def eval(self, sample):
        import RNA
        return RNA.energy_of_struct(sample, self.structure)
    ## @return "EX" where X is the index
    def idstring(self):
        return "".join(map(str,self.identifier))

## @brief Sampler for RedPrint
class RedprintSampler(ir.MultiDimensionalBoltzmannSampler):
    ## @brief Constructor
    ## @param seqlen sequence length
    ## @param structure_strings structures as dot-bracket strings
    ## @param features the features as list or dictionary
    ## @param method method for tree decomposition (0: libhtd, 1-5: respective /strategy/ of TDlib)
    ## @param model base pair or stacking model, sepcified as respective string "bp"/"basepair" or "stack"/"stacking"
    ## @param no_red_constrs no redundant constraints flag (mainly for experimenting)
    def __init__(self, seqlen, structure_strings, features, **kwargs):
        super().__init__( features )

        self.seqlen = seqlen
        self.structures = list(map(rna.parseRNAStructureBps,structure_strings))

        def optarg(feat, default):
            if feat in kwargs:
                return kwargs[feat]
            else:
                return default

        self.method = optarg("method", 0)
        self.model = optarg("model", "bp")
        self.no_red_constrs = optarg("no_red_constrs", False)

        self.setup_engine()

    ## @brief Generate feature record
    ## @param structure_strings structures as dot bracket strings
    ## @param energy_weights
    ## @param gcweight
    ## @param energy_targets
    ## @param gctarget
    ## @param energy_tolerances
    ## @param gctolerance
    @staticmethod
    def gen_features( structure_strings,
                      energy_weights, gcweight,
                      energy_targets=None, gctarget=None,
                      energy_tolerances=None, gctolerance=None ):
        if energy_targets == None:
            energy_targets = [None] * len(energy_weights)

        if energy_tolerances == None:
            energy_tolerances = [None] * len(energy_weights)

        features = [ EnergyFeature(i,
                                   structure_strings[i],
                                   energy_weights[i],
                                   energy_targets[i],
                                   energy_tolerances[i])
                     for i in range(len(energy_weights)) ]
        features.append( GCFeature(gcweight, gctarget, gctolerance) )
        return { f.identifier:f for f in features }

    ## @brief Generate constraint network
    ## @param features dictionary of features
    ## @return constraint network
    def gen_constraint_network(self, features):
        ## build constraint network
        if self.model in ["bp", "basepair"]:
            cn = RNAConstraintNetworkBasePair( self.seqlen,
                                               self.structures,
                                               self.features)
        elif self.model in ["stack", "stacking"]:
            cn = RNAConstraintNetworkStacking( self.seqlen,
                                               self.structures,
                                               self.features )
        else:
            print("Model", self.model,
                  "unknown! Please see help for supported modules.")
            exit(-1)

        return cn

    ## @brief Generate tree decomposition
    ## @param cn constraint network
    ## @return tree decomposition
    def gen_tree_decomposition(self, cn):
        ## make tree decomposition
        return RNATreeDecomposition( cn,
                                     add_red_constrs = not self.no_red_constrs,
                                     method=self.method )

    ## @brief Calculate sample
    ## @return sampled RNA sequence
    def sample(self):
        return ir.values_to_sequence(super().sample().values())


# END Redprint Library
# ##########################################################

# ##########################################################
# main
#

## @brief command line tool definition
## @param args command line arguments
def main(args):
    import RNA

    ## init seed
    if args.seed == None:
        ir.seed(random.randint(0,2**31))
    else:
        ir.seed(args.seed)

    ir.set_bpenergy_table( ir.params_bp )
    ir.set_stacking_energy_table( ir.params_stacking )

    ## read instance
    with open(args.infile) as infh:
        structures = rna.read_inp(infh)

    weights=args.weight
    if weights is None: weights=[1]
    if len(weights) < len(structures):
        weights.extend([weights[-1]]*(len(structures)-len(weights)))

    targets=args.target
    if targets is None: targets=[0]
    if len(targets) < len(structures):
        targets.extend([targets[-1]]*(len(structures)-len(targets)))

    tolerances=args.tolerance
    if tolerances is None: tolerances=[1]
    if len(tolerances) < len(structures):
        tolerances.extend([tolerances[-1]]*(len(structures)-len(tolerances)))


    if len(structures) == 0:
        print("At least one structure is required in input")
        exit(-1)

    seqlen = len(structures[0])

    features = RedprintSampler.gen_features(structures, weights, args.gcweight,
                                            targets, args.gctarget,
                                            tolerances, args.gctolerance)

    sampler = RedprintSampler( seqlen, structures, features,
                               method = args.method,
                               model = args.model,
                               no_red_constrs = args.no_red_constrs )

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
    for seq in sample_generator:

        print(seq,end='')
        if args.turner:
            for i,struc in enumerate(structures):
                feat_id, value = fstats.record( sampler.features[("E",i)], seq )
                print(" {}={:3.2f}".format(feat_id,value),end='')

        if args.gc:
            feat_id, value = fstats.record( sampler.features["GC"], seq )
            print(" GC={:3.2f}".format(value),end='')

        if args.checkvalid:
            for i,struc in enumerate(structures):
                if not rna.is_valid(seq, struc):
                    print(" INVALID{}: {}".format(i,str(rna.invalid_bps(seq,struc))),end='')
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
    parser = argparse.ArgumentParser(description="Boltzmann sampling for RNA design with multiple target structures")
    parser.add_argument('infile', help="Input file")
    parser.add_argument('--method', type=int, default=0,
                        help="Method for tree decomposition (0: use htd; otherwise pass to TDlib as strategy)")
    parser.add_argument('-n','--number', type=int, default=10, help="Number of samples")
    parser.add_argument('--seed', type=int, default=None,
                        help="Seed infrared's random number generator (def=auto)")

    parser.add_argument('-v','--verbose', action="store_true", help="Verbose")
    parser.add_argument('--model', type=str, default="bp",
                        help="Energy model used for sampling [bp=base pair model,stack=stacking model]")
    parser.add_argument('--turner', action="store_true",
                        help="Report Turner energies of the single structures for each sample")
    parser.add_argument('--gc', action="store_true",
                        help="Report GC contents of the single structures for each sample")

    parser.add_argument('--checkvalid', action="store_true",
                        help="Check base pair complementarity for each structure and for each sample")

    parser.add_argument('--no_red_constrs', action="store_true",
                        help="Do not add redundant constraints")

    parser.add_argument('--gcweight', type=float, default=1, help="GC weight")
    parser.add_argument('--weight', type=float, action="append",
                        help="Energy weight (def=1; in case, use last weight for remaining structures)")

    parser.add_argument('--mdbs', action="store_true",
                        help="Perform multi-dim Boltzmann sampling to aim at targets")

    parser.add_argument('--gctarget', type=float, default=50, help="GC target")
    parser.add_argument('--target', type=float, action="append",
                        help="Energy target (def=0; in case, use last target for remaining structures)")

    parser.add_argument('--gctolerance', type=float, default=5, help="GC tolerance")
    parser.add_argument('--tolerance', type=float, action="append",
                        help="Energy tolerance (def=1; in case, use last tolerance for remaining structures)")


    parser.add_argument('--plot_td', action="store_true",
                        help="Plot tree decomposition")

    args=parser.parse_args()

    main(args)
