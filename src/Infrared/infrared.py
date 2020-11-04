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
from libinfrared import seed,Constraint,Function

from treedecomp import TreeDecomposition, TreeDecompositionFactory

import abc
import rna_support as rna


## @brief exception to signal inconsistency, when consistency would be required
class ConsistencyError(RuntimeError):
    def __init__(self, arg):
        self.args = [arg]

## @brief Constraint network base class
##
## The constraint network typically holds the problem instance-specific
## constraints and functions. The fields and methods of this class are
## up to the user; typically, one could define fields dependencies,
## functions, and constraints.
class ConstraintNetwork:
    def __init__(self, *, domains=None, varnum=None, constraints=[], functions=[]):
        if type(domains) is list:
            assert(varnum is None)
            self._domains = domains
        else:
            assert(type(domains) is int)
            self._domains = [domains] * varnum

        self._constraints = constraints
        self._functions = functions

    ## @brief infer dependencies from the functions and constraints
    ##
    ## @param non_redundant whether the dependency list is made non-redundant
    ## @returns list of lists of indices of variables that depend on each other either through functions or constraints
    def get_dependencies(self, non_redundant=True):
        deps = [ x.vars() for x in self.get_functions() + self.get_constraints() ]

        if non_redundant:
            deps = self._remove_redundant_dependencies(deps)

        return deps

    ## @brief list of all functions
    def get_functions(self):
        return self._functions

    ## @brief list of all constraints
    def get_constraints(self):
        return self._constraints

    def get_varnum(self):
        return len(self._domains)

    # list of the domain sizes of each variable
    def get_domains(self):
        return self._domains

    ## @brief Removes redundant dependencies
    ##
    ## @param deps list of dependencies (where a dependency is a list of indices)
    ##
    ## Removes all dependencies that are subsumed by other dependencies in deps
    ##
    ## @return pruned list of dependencies
    @staticmethod
    def _remove_redundant_dependencies(deps):
        def sublist(xs,ys):
            return all(x in ys for x in xs)
        def subsumed(dep,deps):
            return any( len(dep)<len(dep2) and sublist(dep,dep2) for dep2 in deps )
        return [ dep for dep in deps if not subsumed(dep,deps) ]


## @brief Cluster tree (wrapping the cluster tree class of the C++ engine)
class ClusterTree:
    def __init__(self, cn, *, td_factory=TreeDecompositionFactory(), td=None):
        if td is None:
            td = td_factory.create(cn.get_varnum(), cn.get_dependencies())

        self.cn = cn

        self._bagsets = list( map(set, td.get_bags()) )


        self.td = td
        self.ct = self.construct_cluster_tree( cn.get_domains(), td )

    ## @brief Construct the cluster tree object of the C++ engine
    ##
    ## domains can either specify a uniform domain size, or a
    ## list of all domain sizes
    def construct_cluster_tree(self, domains, td):
        bagconstraints, bagfunctions = self.get_bag_assignments()

        ct = libir.ClusterTree(domains);

        # keep record of all non-root nodes
        children = set()

        for bagidx in td.toposorted_bag_indices():
            if not bagidx in children:
                # enumerate subtree
                stack = [(None,bagidx)]
                while stack:
                    (p,i) = stack[-1]
                    stack = stack[:-1]
                    bagvars = sorted(list(self._bagsets[i]))

                    if p==None:
                        cluster = ct.add_root_cluster(bagvars)
                    else:
                        cluster = ct.add_child_cluster(p,bagvars)

                    for x in bagconstraints[i]:
                        ct.add_constraint(cluster, x)

                    for x in bagfunctions[i]:
                        ct.add_function(cluster, x)

                    for j in td.adj[i]:
                        children.add(j)
                        stack.append((cluster,j))
        return ct

    ## @brief evaluates the cluster tree
    ## @return partition function
    ##
    ## @note Evaluation is a potentially (depending on the treewidth) costly operation.
    ## The method does not re-evaluate the tree if this was already done
    def evaluate(self):
        return self.ct.evaluate()

    ## @brief evaluates the cluster tree (and thereby checks consistency)
    ## @return whether the constraints are consistent
    ##
    ## @note does not re-evaluate the tree if this was already done
    def is_consistent(self):
        return self.ct.is_consistent()

    ## @brief generate sample
    ## @returns a raw sample
    ##
    ## @note raises exception ConsistencyError if the network is inconsistent.
    ## If the cluster tree was not evaluated (or consistency checked) before,
    ## it will be evaluated once on-demand.
    def sample(self):
        if not self.is_consistent():
            raise ConsistencyError("Inconsistent constraint model")
        return self.ct.sample()

    def get_td(self):
        return self.td

    def get_bagsets(self):
        return self._bagsets

    ## @brief Get assignments of functions and constraints to the bags
    ##
    ## straightforward non-redundant assignment of constraints and functions,
    ## each to some bag that contains all of their variables
    ##
    ## assumes constraints and functions specified in self.cn
    def get_bag_assignments(self):
        bagconstraints = self.assign_to_bags(self.cn.get_constraints())
        bagfunctions = self.assign_to_bags(self.cn.get_functions())
        return (bagconstraints, bagfunctions)


    ## @brief Get the indices of all bags that contain a set of variables
    ## @param bvars the set of variables
    ## @return list of indices of the bags that contain bvars
    def find_all_bags(self, bvars):
        return [ i for i,bag in enumerate(self._bagsets) if all( x in bag for x in bvars ) ]

    ## @brief Find a bag that contains a set of variables
    ## @param bvars the set of variables
    ## @return index of first bag that contains bvars (or None if there is none)
    def find_bag(self, bvars):
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
        bagconstraints = { i:[]  for i in range(len(self._bagsets)) }
        for cons in constraints:
            bagconstraints[self.find_bag(cons.vars())].append(cons)
        return bagconstraints

    ## @brief assign constraints or functions to all possible bags
    ## @param constraints list of constraints
    ## @return list where constraints are placed at corresponding bag indices
    ##
    ## Assigns constraint/function is to each bags that contains its
    ## dependencies.
    def assign_to_all_bags(self, constraints):
        bagconstraints = { i:[]  for i in range(len(self._bagsets)) }
        for cons in constraints:
            for bidx in self.find_all_bags(cons.vars()):
                bagconstraints[bidx].append(cons)
        return bagconstraints


## @brief Feature in multi-dimensional Boltzmann sampling
##
## A feature is one component of the sampling energy functions, which
## can be targeted due to a dedicated weight. It defines a method
## value, which determines the feature's value for a sample.
##
## A feature defines weight, target value, and tolerance. The latter
## two are used only in case of multi-dimensional Boltzmann sampling,
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
## deviation (or the entire distribution) of each recorded feature.
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
        return ' '.join(["{}={:3.2f} +/-{:3.2f}".format(self.idstring[fid],means[fid],stds[fid]) for fid in self.count])

## @brief Boltzmann sampler (abstract base class)
class BoltzmannSampler:

    ## @brief Construct with features
    ## @param features list or dictionary of the features
    def __init__( self, features, td_factory = TreeDecompositionFactory() ):
        if type(features) == list:
            self.features = { f.identifier:f for f in features }
        else:
            self.features = features

        self._td_factory = td_factory

    ## @brief flag that engine requires setup
    def requires_reinitialization(self):
        self.cn = None
        self.td = None
        self.ct = None

    ## @brief Sets up the constraint network / cluster tree sampling
    ## engine
    ##
    ## @note do nothing, if the engine was already initialized before
    ##
    ## @param skip_ct skip the potentially expensive construction and evaluation of the cluster tree
    def setup_engine(self, *, skip_ct=False):
        # immediately return if ct exists and is not None
        try:
            assert( self.ct != None )
        except:
            pass
        else:
            return

        self.cn = self.gen_constraint_network(self.features)
        self.td = self._td_factory.create(self.cn.get_varnum(), self.cn.get_dependencies())
        if not skip_ct:
            self.ct = self.gen_cluster_tree()

    ## @brief Get the features
    ## @return features
    def get_features(self):
        return self.features

    ## @brief is the network consistent?
    def is_consistent(self):
        self.setup_engine()
        return self.ct.is_consistent()

    ## @brief Plot the tree decomposition to pdf file
    def plot_td(self, dotfilename):
        self.setup_engine(skip_ct = True)
        with open(dotfilename,"w") as dot:
            self.td.writeTD(dot)
        import treedecomp
        treedecomp.dotfile_to_pdf(dotfilename)
        os.remove(dotfilename)

    ## @brief Get tree width
    ## @return tree width
    def treewidth(self):
        self.setup_engine(skip_ct = True)
        return self.td.treewidth()

    ## @brief Compute next sample
    ## @return sample
    def sample(self):
        self.setup_engine()
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

    ## @brief Generate the populated cluster tree
    ## @param td tree decomposition
    ## @return cluster tree
    def gen_cluster_tree(self):
        ## make cluster tree
        return ClusterTree(self.cn, td = self.td)

## @brief Multi-dimensional Boltzmann sampler (abstract base class)
class MultiDimensionalBoltzmannSampler(BoltzmannSampler):
    def __init__( self, features, td_factory=TreeDecompositionFactory() ):
        super().__init__( features, td_factory )

        self.samples_per_round = 200
        self.tweak_base = 1.01

    ## @brief whether the sample is of good quality
    ## @param features dictionary of features
    ##
    ## checks whether the sample approximately meets the targets
    def is_good_sample(self, features):
        for f in features.values():
            if abs(f.value - f.target) > f.tolerance:
                return False
        return True

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

        counter = 0
        while True:
            fstats = FeatureStatistics()
            for i in range(self.samples_per_round):
                counter+=1
                sample = self.sample()
                returned_features = fstats.record( self.features, sample )

                if self.is_good_sample(returned_features):
                    # print(counter)
                    yield sample

            last_means=means
            means = fstats.means()

            # modify weight of each feature
            for fid,f in self.features.items():
                f.weight = f.weight *  self.tweak_base**(means[fid] -f.target)

            self.requires_reinitialization()

