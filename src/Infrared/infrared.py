#!/usr/bin/env python3

# -----------------------------
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

import itertools
import os

import libinfrared as libir
from libinfrared import *

import treedecomp
import abc
import rna_support as rna


## @brief Constraint network base class
##
## The consraint network typically holds the problem instance-specific
## constraints and functions. The fields and methods of this class are
## up to the user; typically, one could define fields dependencies,
## functions, and constraints.
class ConstraintNetwork:
    def __init__(self):
        pass

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

## @brief Base class for tree decomposition defining some more
## commonly useful methods.
##
## To see an example for the use of this class, check out the redprint
## tool.
##
class TreeDecomposition:
    ## @brief Constructor
    ## @param varnum number of variables
    ## @param dependencies list of dependencies
    ## @param method tree decomposition method
    def __init__(self, varnum, dependencies, *, method=0):
        self.varnum = varnum

        # from dependencies generate list of binary edges
        bindependencies  = self.expand_to_cliques(dependencies)

        # generate tree decomposition -> bags, edges
        self.td = treedecomp.makeTD(varnum, bindependencies, method = method)

        self.bagsets = list(map(set,self.td.bags))

        self.domains = None

    ## @brief Get assignments of functions and constraints to the bags
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

    ## @brief Construct the cluster tree
    ##
    ## Requires deriving classes to specialize self.domains and
    ## self.get_bag_assignments()
    ##
    ## self.domains can either specify a uniform domain size, or a
    ## list of all domain sizes
    def construct_cluster_tree(self):

        bagconstraints, bagfunctions = self.get_bag_assignments()

        if type(self.domains) == int:
            ct = libir.ClusterTree(self.varnum, self.domains);
        else:
            ct = libir.ClusterTree(self.domains);

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


## @brief Feature in multi-dimensional Boltzmann sampling
##
## A feature is one component of the sampling energy functions, which
## can be targeted due to a dedicated weight. It defines a method
## value, which determines the feature's value for a sample.
##
## A feature defines weight, target value, and tolerance. The latter
## two are used only in case of multi-dimenstional Boltzmann sampling,
## which modifies the features weight based on the difference between
## the estimated mean feature value and the target value.
##
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


## @brief Keeping statistics on features
##
## This class allows recording values of multiple features for a
## series of samples; it can be queried for mean and standard
## deviation (or the entire distribution) of each recored feature.
class FeatureStatistics:
    ## @brief constructor
    ## @param keep Keep all recorded features in memory
    def __init__(self, keep=False):
        self.keep = keep

        self.idstring = dict()
        self.count = dict()
        self.sums = dict()
        self.sqsums = dict()
        self.features = dict()

    ## @brief Record feature values
    ## @param feature a single feature or a dictionary of features
    ## @param sample a sample that can be evaluated by the feature(s)
    ## @return pair of feature id string and value
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

    ## @brief check whether any features have been recorded
    ## @return whether empty (since no feature has been recorded)
    def empty(self):
        return not self.count

    ## @brief means of recorded features
    ## @return dictionary of means  (at feature ids as keys)
    def means(self):
        return { k : self.sums[k]/self.count[k] for k in self.sums }

    ## @brief variances of recorded features
    ## @return dictionary of variances (at feature ids as keys)
    def variances(self):
        means = self.means()
        return { k : self.sqsums[k]/self.count[k] - means[k]**2 for k in self.sqsums }

    ## @brief standard deviations of recorded features
    ## @return dictionary of standard deviations (at feature ids as keys)
    def stds(self):
        return { k:val**0.5 for k,val in self.variances().items() }

    ## @brief Report features to standard output
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
    ## @brief Construct with features
    ## @param features list or dictionary of the features
    def __init__(self, features):
        if type(features) == list:
            self.features = { f.identifier:f for f in features }
        else:
            self.features = features

    ## @brief Sets up the constraint network / cluster tree sampling
    ## engine
    def setup_engine(self):
        self.cn = self.gen_constraint_network(self.features)
        self.td = self.gen_tree_decomposition(self.cn)
        self.ct = self.gen_cluster_tree(self.td)
        self.requires_evaluation = True

    ## @brief Get the features
    ## @return features
    def get_features(self):
        return self.features

    ## @brief Plot the tree decomposition to pdf file
    def plot_td(self, dotfilename):
        with open(dotfilename,"w") as dot:
            self.td.td.writeTD(dot)
        treedecomp.dotfile_to_pdf(dotfilename)
        os.remove(dotfilename)

    ## @brief Get tree width
    ## @return tree width
    def treewidth(self):
        return self.td.td.treewidth()

    ## @brief Compute next sample
    ## @return sample
    def sample(self):
        if self.requires_evaluation:
            self.ct.evaluate()
            self.requires_evaluation = False

        return self.ct.sample()

    ## @brief Sample generator
    def samples(self):
        while(True):
            yield self.sample()

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

        self.samples_per_round = 200
        self.tweak_base = 1.01

        self.counter_total_samples = 0
        self.counter_good_samples = 0

    ## @brief whether the sample is of good quality
    ## @param features dictionary of features
    ##
    ## checks whether the sample approximately meets the targets
    def is_good_sample(self, features):
        for f in features.values():
            if abs(f.value - f.target) > f.tolerance:
                return False
        return True

    def success_rate(self):
        return self.counter_good_samples / float(self.counter_total_samples)

    def reset_sample_counter(self):
        self.counter_good_samples = 0
        self.counter_total_samples = 0

    ## @brief Generator of targeted samples
    ##
    ## Performs multi-dimensional Boltzmann sampling: every
    ## self.samples_per_round many samples, the feature means are
    ## estimated and the weights are recalibrated. Each generated
    ## sample is tested for falling into target +/- tolerance for all
    ## features, in which case it is yielded.
    ##
    ## self.tweak_base controls the speed of recalibration of weights due to
    ## the formula
    ## ```
    ## weight = weight * tweak_base**(mean -target)
    ## ```
    def targeted_samples(self):
        means=None

        while True:
            self.setup_engine()
            fstats = FeatureStatistics()
            for i in range(self.samples_per_round):
                self.counter_total_samples += 1
                sample = self.sample()
                returned_features = fstats.record( self.features, sample )

                if self.is_good_sample(returned_features):
                    self.counter_good_samples += 1
                    yield sample

            last_means=means
            means = fstats.means()

            # modify weight of each feature
            for fid,f in self.features.items():
                f.weight = f.weight *  self.tweak_base**(means[fid] -f.target)

# ------------------------------------------------------------
# RNA specific definitions

## @brief Convert integer (variable value) to nucleotide
## @note encoding A=0, C=1, G=2, U=3
def value_to_nucleotide(x):
    return "ACGU"[x]

## @brief Convert list of integers (variable values) to string
## (sequence) of nucleotides
def values_to_sequence(xs):
    return "".join(map(value_to_nucleotide, xs))

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

## @brief set the bp energy table for Infrared
##
## @param params dictionary of parameters or a table of the parameters
## as expected by infrared::rnadesign
def set_bpenergy_table(params):
    if type(params) == dict:
        params = list(map(lambda x: params[x],
                          [ "AU_IN", "GC_IN", "GU_IN",
                            "AU_TERM", "GC_TERM", "GU_TERM" ] ))
    libir.BPEnergy.set_energy_table(params)

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
    libir.StackEnergy.set_energy_table(params)
