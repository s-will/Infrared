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

import abc

## Parameters for the base pair model (magic params from the Redprint paper)
params_bp = { "GC_IN": -2.10208, "AU_IN": -0.52309, "GU_IN": -0.88474,
              "GC_TERM": -0.09070, "AU_TERM": 1.26630, "GU_TERM": 0.78566 }

## Parameters for the stacking model (magic params from the Redprint paper)
params_stacking = { "AUAU": -0.18826, "AUCG": -1.13291, "AUGC": -1.09787,
                    "AUGU": -0.38606, "AUUA": -0.26510, "AUUG": -0.62086,
                    "CGAU": -1.11752, "CGCG": -2.23740, "CGGC": -1.89434,
                    "CGGU": -1.22942, "CGUA": -1.10548, "CGUG": -1.44085,
                    "GUAU": -0.55066, "GUCG": -1.26209, "GUGC": -1.58478,
                    "GUGU": -0.72185, "GUUA": -0.49625, "GUUG": -0.68876 }

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
##
## @param params dictionary of parameters or a table of the parameters
## as expected by infrared::rnadesign
def set_bpenergy_table(params):
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          [ "AU_IN", "GC_IN", "GU_IN",
                            "AU_TERM", "GC_TERM", "GU_TERM" ] ))
    ir.BPEnergy.set_energy_table(params)

## @brief set the stacking energy table for Infrared
##
## @param params dictionary of parameters or a table of the parameters
## as expected by infrared::rnadesign
def set_stacking_energy_table(params):
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          [ "AUAU", "AUUA",
                            "AUCG", "AUGC",
                            "AUGU", "AUUG",

                            "CGAU", "CGUA",
                            "CGCG", "CGGC",
                            "CGGU", "CGUG",

                            "GUAU", "GUUA",
                            "GUCG", "GUGC",
                            "GUGU", "GUUG" ] ))
    ir.StackEnergy.set_energy_table(params)

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
        gc_control_funs = [ ( [i], [ir.GCControl(i,self.gcweight)] )
                            for i in range(self.seqlen) ]
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

## @brief Construct and hold constraint network for mult-target design
## based in the base pair energy model
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
                                       [ir.BPEnergy(i,j, not (i-1,j+1) in structureset,
                                                    self.weights[k])] )
                                     for (i,j) in structure ] )

        self.add_gc_control()

## @brief Construct and hold constraint network for mult-target design
## based in the stacking energy model
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

## @brief Base class for tree decomposition defining some more
## commonly useful methods
##
class TreeDecomposition:
    def __init__(self, varnum, dependencies, *, method=0):
        self.varnum = varnum

        # from dependencies generate list of binary edges
        bindependencies  = self.expand_to_cliques(dependencies)

        # generate tree decomposition -> bags, edges
        self.td = treedecomp.makeTD(varnum, bindependencies, method = method)

        self.bagsets = list(map(set,self.td.bags))

        self.domains = None

    @abc.abstractmethod
    def get_bag_assignments(self):
        pass

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

    ## @brief assign functions and constraints to bags
    ##
    ## @return bagconstraints, bagfuntions
    ##
    ## the function is called by construct_cluster_tree, must be overriden
    @abc.abstractmethod
    def make_bag_assignments(self):
        pass

    ## @brief Construct the cluster tree
    ##
    ## Requires deriving classes to specialize self.domains and
    ## self.make_bag_assignments()
    ##
    ## self.domains can either specify a uniform domain size, or a
    ## list of all domain sizes
    def construct_cluster_tree(self):

        bagconstraints, bagfunctions = self.get_bag_assignments()

        if type(self.domains) == int:
            ct = ir.ClusterTree(self.varnum, self.domains);
        else:
            ct = ir.ClusterTree(self.domains);

        # keep record of all non-root nodes
        children = set()

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

                    for x in bagconstraints[i]:
                        ct.add_constraint(cluster, x)

                    for x in bagfunctions[i]:
                        ct.add_function(cluster, x)

                    for j in self.td.adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct


## @brief Holds tree decomposition for RNA design, constructs RNA design cluster tree
##
class RNATreeDecomposition(TreeDecomposition):

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

## @brief Feature in multi-dimensional Boltzmann sampling
##
## A feature is one component of the sampling energy functions,
## which can be targeted due to a dedicated weight.
class Feature:
    def __init__(self, identifier, weight, target=None, tolerance=None):
        self.identifier = identifier
        self.weight = weight
        self.target = target
        self.tolerance = tolerance

    ## @brief Printable identifier
    def idstring(self):
        return str(self.identifier)

    ## @brief Evaluate feature for given sample
    @abc.abstractmethod
    def eval(self, sample):
        return
## @brief GC content feature
class GCFeature(Feature):
    def __init__(self, weight, target=None, tolerance=None):
        super().__init__( "GC", weight, target, tolerance)
    def eval(self, sample):
        return rna.GC_content(sample) * 100

## @brief Turner energy feature
class EnergyFeature(Feature):
    def __init__(self, index, structure, weight, target=None, tolerance=None):
        super().__init__( ("E",index), weight, target, tolerance )
        self.structure = structure
    def eval(self, sample):
        import RNA
        return RNA.energy_of_struct(sample, self.structure)
    def idstring(self):
        return "".join(map(str,self.identifier))

## @brief Keeping statistics on features
##
## @param keep Keep all recorded features in memory
class FeatureStatistics:
    def __init__(self, keep=False):
        self.keep = keep

        self.idstring = dict()
        self.count = dict()
        self.sums = dict()
        self.sqsums = dict()
        self.features = dict()

    def record(self, feature, sample):
        
        if type(feature) == dict:
            for k in feature:
                fid, value = self.record(feature[k], sample)
                feature[k].value = value
            return feature
        
        fid = feature.identifier
        value   = feature.eval(sample)

        if fid not in self.count:
            self.idstring[fid] = feature.idstring()
            self.count[fid] = 0
            self.sums[fid] = 0
            self.sqsums[fid] = 0
            if self.keep:
                self.features[fid] = []

        self.count[fid] += 1
        self.sums[fid] += value
        self.sqsums[fid] += value**2

        if self.keep:
            self.features[fid].append(value)

        return feature.idstring(), value

    def empty(self):
        return not self.count

    ## @brief means of recorded features
    def means(self):
        return { k : self.sums[k]/self.count[k] for k in self.sums }

    ## @brief variances of recorded features
    def variances(self):
        means = self.means()
        return { k : self.sqsums[k]/self.count[k] - means[k]**2 for k in self.sqsums }

    ## @brief standard deviations of recorded features
    def stds(self):
        return { k:val**0.5 for k,val in self.variances().items() }

    def report(self):
        means = self.means()
        stds= self.stds()
        for fid in self.count:
            print(" {}={:3.2f} +/-{:3.2f}"
                  .format(self.idstring[fid],means[fid],stds[fid]),end='')
        print()

## @brief Boltzmann sampler (abstract base class)
##
## derived classes must override gen_cluster_tree
class BoltzmannSampler:
    def __init__(self, features):
        if type(features) == list:
            self.features = { f.identifier:f for f in features }
        else:
            self.features = features

    def setup_engine(self):
        self.cn = self.gen_constraint_network(self.features)
        self.td = self.gen_tree_decomposition(self.cn)
        self.ct = self.gen_cluster_tree(self.td)
        self.requires_evaluation = True

    def get_features(self):
        return self.features

    def plot_td(self, dotfilename):
        with open(dotfilename,"w") as dot:
            self.td.td.writeTD(dot)
        treedecomp.dotfile_to_pdf(dotfilename)
        os.remove(dotfilename)

    def treewidth(self):
        return self.td.td.treewidth()

    def sample(self):
        if self.requires_evaluation:
            self.ct.evaluate()
            self.requires_evaluation = False

        return self.ct.sample()

    ## @brief Generate the constraint network
    ## @param features the features (containing their weights)
    ## @return the constraint network
    @abc.abstractmethod
    def gen_constraint_network(self, features):
        pass

    ## @brief Generate the tree decomposition
    ## @param cn the constraint network
    ## @return the tree decomposition
    @abc.abstractmethod
    def gen_tree_decomposition(self, cn):
        pass

    ## @brief Generate the populated cluster tree
    ## @param td tree decomposition
    ## @return cluster tree
    def gen_cluster_tree(self, td):
        ## make cluster tree
        ct = td.construct_cluster_tree()
        return ct

## @brief Multi-dimenstional Boltzmann sampler (abstract base class)
class MultiDimensionalBoltzmannSampler(BoltzmannSampler):
    def __init__(self, features):
        super().__init__(features)
        
        self.samples_per_round = 100
        self.tweak_base = 1.01
        
    ## @brief whether the sample is of good quality
    ##
    ## checks whether the sample approximately meets the targets
    def is_good_sample(self, features):
        for f in features.values():
            if abs(f.value - f.target) > f.tolerance:
                return False
        return True
        
    def targeted_samples(self):
        means=None
        while True:
            self.setup_engine()
            fstats = FeatureStatistics()
            for i in range(self.samples_per_round):
                sample = self.sample()
                returned_features = fstats.record( self.features, sample )
                
                if self.is_good_sample(returned_features):
                    yield sample
                
            last_means=means
            means = fstats.means()
            
            # modify weight of each feature
            for fid,f in self.features.items():
                f.weight = f.weight *  self.tweak_base**(means[fid] -f.target)
                # print(" {} = {:3.2f}->{:3.2f} ({:3.2f})".format(f.idstring(), means[fid], f.target, f.weight))
               
            # print("==============================")

## @brief Sampler for RedPrint
class RedprintSampler(MultiDimensionalBoltzmannSampler):
    def __init__(self, seqlen, structure_strings, features, **kwargs):
        super().__init__( features )

        self.seqlen = seqlen
        self.structures = list(map(rna.parseRNAStructureBps,structure_strings))
        
        def optarg(feat, default):
            if hasattr(kwargs,feat):
                return kwargs[feat]
            else:
                return default

        self.method = optarg("method", 0)
        self.model = optarg("model", "bp")
        self.no_red_constrs = optarg("no_red_constrs", False)

        self.setup_engine()

    ## @brief Generate feature record
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
    def gen_constraint_network(self, features):
        weights = [ features[("E",i)].weight for i in range(len(self.structures)) ]
        gcweight = features["GC"].weight

        ## build constraint network
        if self.model in ["bp", "basepair"]:
            cn = RNAConstraintNetworkBasePair( self.seqlen,
                                               self.structures,
                                               weights,
                                               gcweight)
        elif self.model in ["stack", "stacking"]:
            cn = RNAConstraintNetworkStacking( self.seqlen,
                                               self.structures,
                                               weights,
                                               gcweight )
        else:
            print("Model", self.model,
                  "unknown! Please see help for supported modules.")
            exit(-1)

        return cn

    def gen_tree_decomposition(self, cn):
        ## make tree decomposition
        return RNATreeDecomposition( cn,
                                     add_red_constrs = not self.no_red_constrs,
                                     method=self.method )

    def sample(self):
        return values2seq(super().sample().values())

    ## @brief sample generator
    def samples(self):
        while(True):
            yield values2seq(super().sample().values())

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

    set_bpenergy_table( params_bp )
    set_stacking_energy_table( params_stacking )

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

    fstats = FeatureStatistics()

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
