#!/usr/bin/env python3

# -----------------------------
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

import os
import abc
import random
import re
import sys
import inspect
import math

from . import libinfrared as libir

from treedecomp import TreeDecompositionFactory, dotfile_to_pdf, dotfile_to_png
from treedecomp import seed as tdseed

from abc import ABC, abstractmethod

def seed(seed):
    """!@brief seed random number generator of libinfrared and treedecomp

    This seeds the RNG of lib infrared (C++ side) and as well
    the random number generator used by randomization in the TreeDecomposition,
    both with the same number

    This does not seed Python's global RNG

    @param seed integer used as seed

    Without argument or seed==None, use pythons built-in random.seed() to generate
    a seed.
    """

    if seed is None:
        random.seed()
        seed = random.randint(0,2**31)

    libir.seed(seed)
    tdseed(seed)


## @brief exception to signal inconsistency, when consistency would be required
class ConsistencyError(RuntimeError):
    def __init__(self, arg):
        self.args = [arg]

class EvaluationAlgebra(ABC):
    """!@brief Algebra for evaluating a constraint network
    """

    @staticmethod
    @abstractmethod
    def cluster_tree(*args,**kwargs):
        """!@brief the infrared cluster tree
        @return the cluster tree that evaluates under the specific algebra"
        """
        pass

    @staticmethod
    @abstractmethod
    def value( weight, value ):
        """!@brief Combine weight and value of weighted function
        @param weight the weight
        @param value the value

        @return combination of weight and value according to concrete algebra
        """
        pass

class PFEvaluationAlgebra(EvaluationAlgebra):
    """!@brief Partition function algebra for sampling"""
    def cluster_tree(*args,**kwargs):
        return libir.PFClusterTree(*args,**kwargs)

    def value( weight, value ):
        return math.exp( weight * value )

class ArcticEvaluationAlgebra(EvaluationAlgebra):
    """!@brief Maximization algebra for optimization"""

    def cluster_tree(*args,**kwargs):
        return libir.ArcticClusterTree(*args,**kwargs)

    def value( weight, value ):
        return weight * value

class WeightedFunction(libir.Function):
    """!@brief function of a constraint network

    WeightedFunction have methods value() and weight(); value depends
    on the variables defined at construction and returned by vars().
    """

    _algebra = PFEvaluationAlgebra

    def __init__( self, variables ):
        super().__init__(variables)
        self._weight = 1

    @abstractmethod
    def value(self):
        pass

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, weight):
        self._weight = weight

    @staticmethod
    @property
    def algebra(algebra):
        _algebra = algebra

    def __call__(self, a):
        return self._algebra.value( self.weight, self.value(a) )


def _generic_def_function_class( classname, init, value, module="__main__",
            parentclass = WeightedFunction, valuefunname = "value" ):
    """!@brief Create a class of infrared weighted functions or constraints
    @param init function to create dependency list from constructor argument(s)
    @param value function to compute the weighted functions's value from assignment values
    @param env environment, in which the class is created (users could e.g. pass env=locals())

    Defines a new class with name 'classname' (by default in the module's namespace)

    Objects of the defined class can be used in the construction of the infrared constraint network.
    Note that `init` defines the dependencies in the order of the (positional) arguments of `value`.

    @note the value function can depend on arguments to the function init,
    which will be automatically stored in the class and passed on.
    """

    def _init(self, *args, **kwargs):

        variables = init(*args,**kwargs)
        super(self.__class__, self).__init__( variables )

        siginit = inspect.signature(init)
        sigvalue = inspect.signature(value)
        for i,kw in zip( range(len(args)), siginit.parameters ):
            kwargs[kw] = args[i]
        self._args = { k:kwargs[k] for k in kwargs if k in sigvalue.parameters }

    def _value(self,a):
        a = a.values()
        params = [ a[var] for var in self.vars() ]
        return value( *params, **self._args )

    def _str(self):
        return '{} on {}'.format(self.__class__, self._vars)

    newclass = type( classname, (parentclass,),
                {
                 "__init__": _init,
                 valuefunname: _value,
                 "__str__": _str,
                 "name": lambda self: classname
                 })

    sys.modules[module].__dict__[classname] = newclass

def def_function_class( classname, init, value, module="__main__" ):
    _generic_def_function_class( classname, init, value, module, WeightedFunction, "value" )

def def_constraint_class( classname, init, value, module="__main__" ):
    _generic_def_function_class( classname, init, value, module, libir.Constraint, "__call__" )

class Model:
    """!@brief A constraint model
    """
    def __init__( self ):
        self._constraints = []
        self._functions = dict()
        self._domains = dict()

        self._features = dict()

    def add_variables( self, number, domain, name = 'X' ):
        """!@brief add variable domains
        @param number number of variables
        @param domain domain size; defines the domain values as 0..domain-1
        @param name assign a name to the variable(s)
        """

        if name not in self._domains:
            self._domains[name] = []

        self._domains[name].extend( [ libir.FiniteDomain(domain) for _ in range(number) ] )

    def restrict_domains( self, vars, domain ):
        """!@brief restrict the domain of a variable
        @param vars variable or list of variables, each specified by (name,index); or simply index, then addressing ('X',index)
        @param domain the domain

        @note the domain bounds must be stricter than the original domain
        """
        newdom = libir.FiniteDomain(domain)
        if type(vars) != list:
            vars = [vars]
        for v in vars:
            name,i = v if type(v) == tuple else ('X',v)

            assert(self._domains[name][i].lb() <= newdom.lb())
            assert(self._domains[name][i].ub() >= newdom.ub())
            self._domains[name][i] = newdom

    def add_constraints(self, constraints):
        """!@brief add constraints to the model
        @param constraints an iterable of constraints or a single constraint
        """
        if hasattr(constraints, '__iter__'):
            constraints = list(constraints)
        else:
            constraints = [constraints]

        self._constraints.extend(constraints)

    def add_functions( self, functions, group = '*' ):
        """!@brief add constraints to the model
        @param constraints an iterable of constraints or a single constraint
        """

        if group not in self._functions:
            self._functions[group] = []

        # reserve auto feature entry for group or invalidate cached feature
        self._features[group] = None

        if hasattr(functions, '__iter__'):
            functions = list(functions)
        else:
            functions = [functions]

        self._functions[group].extend(functions)

    def num_variables( self, name='X' ):
        """@brief number of variables
        @param name name of the variables series
        @return number of variables of the given series
        """
        if name not in self._domains: return 0
        return len( self._domains[name] )

    def num_all_variables( self ):
        """@brief number of all variables
        @return number of all variables
        """
        return sum( len( self._domains[name] ) for name in self._domains )

    @property
    def domains( self ):
        """
        @brief Domains of the model

        @return A list of the domainsizes of all variables
        """
        return self._domains

    @property
    def constraints( self ):
        """
        @brief Functions of the model
        @return specification of the model's functions
        """
        return self._constraints

    @property
    def functions( self ):
        """
        @brief Functions of the model
        @return specification of the model's functions

        Functions can be either specified as list or dictionary of function groups.
        If defined as list, all functions will be assigned to a single group "Total".
        Functions can be assigned to more than one group.
        @todo: remember to set weights and make functions unique before feeding them to the CN

        Function groups are used to derive the features of the model.
        """
        return self._functions

    @property
    def all_functions( self ):
        """!@brief All functions of the model
        @return list of all functions
        """
        fs = []
        for k in self._functions:
            fs.extend( self._functions[k] )
        return fs

    @property
    def all_domains( self ):
        """!@brief list of all domain indices
        """
        doms = []
        for k in self._domains:
            doms.extend( self._domains[k] )
        return doms

    def _automatic_feature( self, group ):
        """!@brief Automatic feature for function group
        @return the feature

        Feature with value that is derived as sum of the feature functions
        """
        def eval_fun(sample):
            return sum( f.value(sample) for f in self.functions[group] )

        return Feature(group, eval_fun, group = group)

    def add_feature( self, name, group, eval_fun ):
        """!@brief Add a (custom) feature
        @param name name of the feature
        @param group one or several groups of feature controlled functions
        @param eval_fun function to evaluate the feature at an assignment
        """
        self._features[name] = Feature(name, eval_fun, group = group)

    @property
    def features( self ):
        """!@brief Features of the model
        @return dictionary of features
        """
        for name in self._features:
            if self._features[name] is None:
                self._features[name] = self._automatic_feature( name )

        return self._features

    def eval_feature( self, assignment, name ):
        """!@brief evaluate named feature at assignment
        @param assignment the assignment
        @param name name of the feature
        @return value of the feature
        """
        return self.features[name].eval( assignment )

    def set_feature_weight( self, weight, name ):
        """!@brief set the weight of a feature and its function group
        @param weight the weight
        @param name the feature name
        """
        self.features[name].weight = weight
        groups = self.features[name].group
        if type(groups) != list:
            groups = [groups]
        for group in groups:
            for f in self._functions[group]:
                f.weight = weight

    def idx( self, variables ):
        """!@brief raw indices of named variables
        @param variables list of (potentially) names variables; default name 'X'
        @return list of (internal) variable indices
        """
        def convert(var):
            try:
                (name,idx) = var
            except:
                return var
            offset = 0
            for k in self._domains:
                if k==name:
                    break
                offset += len(self._domains[k])
            return offset + idx

        variables = [ convert(var) for var in variables ]
        return variables

## @brief Constraint network
## @todo: merge with Model
##
## The constraint network typically holds the problem instance-specific
## constraints and functions. The fields and methods of this class are
## up to the user; typically, one could define fields dependencies,
## functions, and constraints.
class ConstraintNetwork:
    def __init__(self, model):
        self._domains = model.all_domains
        self._constraints = model.constraints
        self._functions = model.all_functions

    ## @brief infer dependencies from the functions and constraints
    ##
    ## @param non_redundant whether the dependency list is made non-redundant
    ## @returns list of lists of indices of variables that depend on each other either through functions or constraints
    def dependencies(self, non_redundant=True):
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
    def __init__(self, cn, *, td_factory = TreeDecompositionFactory(), td = None, EvalAlg = PFEvaluationAlgebra):
        if td is None:
            td = td_factory.create(cn.get_varnum(), cn.dependencies())

        self._EvalAlg = EvalAlg

        self._cn = cn

        self._bagsets = list( map(set, td.get_bags()) )


        self._td = td
        self._ct = self.construct_cluster_tree( cn.get_domains(), td )

    @property
    def cn(self):
        return self._cn

    @property
    def td(self):
        return self._td

    @property
    def ct(self):
        return self._ct

    ## @brief Construct the cluster tree object of the C++ engine
    ##
    ## @param domains description of the domains 
    def construct_cluster_tree(self, domains, td):
        bagconstraints, bagfunctions = self.get_bag_assignments()

        ct = self._EvalAlg.cluster_tree(domains);

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
        return self._ct.evaluate()

    ## @brief evaluates the cluster tree (and thereby checks consistency)
    ## @return whether the constraints are consistent
    ##
    ## @note does not re-evaluate the tree if this was already done
    def is_consistent(self):
        return self._ct.is_consistent()

    ## @brief generate sample
    ## @returns a raw sample
    ##
    ## @note raises exception ConsistencyError if the network is inconsistent.
    ## If the cluster tree was not evaluated (or consistency checked) before,
    ## it will be evaluated once on-demand.
    def sample(self):
        if not self.is_consistent():
            raise ConsistencyError("Inconsistent constraint model")
        return self._ct.sample()

    def get_td(self):
        return self._td

    def get_bagsets(self):
        return self._bagsets

    ## @brief Get assignments of functions and constraints to the bags
    ##
    ## straightforward non-redundant assignment of constraints and functions,
    ## each to some bag that contains all of their variables
    ##
    ## assumes constraints and functions specified in self._cn
    def get_bag_assignments(self):
        bagconstraints = self.assign_to_bags(self._cn.get_constraints())
        bagfunctions = self.assign_to_bags(self._cn.get_functions())
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


class Feature:
    """!@brief Feature in multi-dimensional Boltzmann sampling
    A feature is one component of the sampling energy functions, which
    can be targeted due to a dedicated weight. It defines a method
    value, which determines the feature's value for a sample.

    A feature defines weight, target value, and tolerance. The latter
    two are used only in case of multi-dimensional Boltzmann sampling,
    which modifies the features weight based on the difference between
    the estimated mean feature value and the target value.
    Moreover, a feature controls one or several groups of functions.
    """
    def __init__(self, identifier, eval_function, weight = 0, target = None, tolerance = None, group = []):
        """@brief Construct feature
        """
        self.identifier = identifier
        self.weight = weight
        self.target = target
        self.tolerance = tolerance
        self.group = group
        self.eval_function = eval_function

    ## @brief Printable identifier
    def idstring(self):
        return str(self.identifier)

    ## @brief Evaluate feature for given sample
    def eval(self, sample):
        return self.eval_function(sample)


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
    ## @param features a dictionary of features
    ## @param values a corresponding dictionary of the feature values
    ## @param sample a sample that can be evaluated by the feature(s)
    ## @return pair of feature id string and value
    def record_features(self, features, values):
        if not type(features) is dict:
            features = { 'dummy': features }
            values = { 'dummy': values }

        for k in features:
            value = values[k]
            fid = features[k].identifier

            if fid not in self.count:
                self.idstring[fid] = features[k].idstring()
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

## @brief Boltzmann sampler
## @todo simplify class, remove methods that were introduced to be overridable (which is not
## needed anymore)
##
## @todo define how precisely BoltzmannSampler is distinguished from ClusterTree !
## (considering that ClusterTree already offers sampling functionality;
## probably, it should not or we rename to PFClusterTree... !
## this gets important, especially when algebras shall be integrated)
##
class BoltzmannSampler:

    ## @brief Construct from model
    ## @param model the constraint model
    def __init__( self, model, td_factory = TreeDecompositionFactory() ):
        self._td_factory = td_factory
        self._model = model
        self.features = model.features

        self.requires_reinitialization()

    ## @brief flag that engine requires setup
    def requires_reinitialization(self):
        self._cn = None
        self._td = None
        self._ct = None

    @property
    def cn(self):
        return self._cn

    @property
    def td(self):
        return self._td

    @property
    def ct(self):
        return self._ct

    @property
    def model(self):
        return self._model

    ## @brief Sets up the constraint network / cluster tree sampling
    ## engine
    ##
    ## @note do nothing, if the engine was already initialized before
    ##
    ## @param skip_ct skip the potentially expensive construction and evaluation of the cluster tree
    def setup_engine(self, *, skip_ct=False):
        # immediately return if ct exists and is not None
        if self._ct != None:
            return

        self._cn = self.gen_constraint_network()

        if self._td is None:
            self._td = self._td_factory.create(self._cn.get_varnum(), self._cn.dependencies())

        if not skip_ct:
            self._ct = self.gen_cluster_tree()

    def is_consistent(self):
        self.setup_engine()
        return self._ct.is_consistent()

    ## @brief Plot the tree decomposition to pdf file
    ## @param filename write to filename
    ## @param to target format, support conversion to "pdf" or "png". Anything else writes graphviz dot format
    def plot_td(self, filename, to="pdf"):
        conversions = {
                "pdf": dotfile_to_pdf,
                "png": dotfile_to_png
                }

        if to in conversions:
            filename = re.sub(f".{to}$","",filename)+".dot"

        self.setup_engine(skip_ct = True)
        with open(filename,"w") as dot:
            self._td.writeTD(dot)

        if to in conversions:
            conversions[to](filename)
            os.remove(filename)

    ## @brief Get tree width
    ## @return tree width
    def treewidth(self):
        self.setup_engine(skip_ct = True)
        return self._td.treewidth()

    ## @brief Compute next sample
    ## @return sample
    def sample(self):
        self.setup_engine()
        return self._ct.sample()

    ## @brief Sample generator
    def samples(self):
        while(True):
            yield self.sample()

    ## @brief Generate the constraint network
    ## @param features the features (containing their weights)
    ## @return the constraint network
    def gen_constraint_network( self ):
        return ConstraintNetwork( self._model )

    ## @brief Generate the populated cluster tree
    ## @param td tree decomposition
    ## @return cluster tree
    def gen_cluster_tree(self):
        ## make cluster tree
        return ClusterTree(self._cn, td = self._td)

## @brief Multi-dimensional Boltzmann sampler
class MultiDimensionalBoltzmannSampler(BoltzmannSampler):
    def __init__( self, model, td_factory=TreeDecompositionFactory() ):
        super().__init__( model, td_factory )

        self.samples_per_round = 1000

        self.tweak_factor = 0.01

    ## @brief whether the sample is of good quality
    ## @param features dictionary of features
    ##
    ## checks whether the sample approximately meets the targets;
    ## check only the targeted features (which have value, target and tolerance)
    def is_good_sample(self, features, values):
        ret = True
        for k,f in features.items():
            try:
                if abs( values[k] - f.target ) > f.tolerance:
                    ret = False
                    break
            except:
                pass
        return ret

    def set_target( self, target, tolerance, featureid ):
        """!@brief Set target of a feature
        @todo fix details
        """
        f = self.features[featureid]
        f.target = target
        f.tolerance = tolerance

    ## @brief Generator of targeted samples
    ##
    ## Performs multi-dimensional Boltzmann sampling: every
    ## self.samples_per_round many samples, the feature means are
    ## estimated and the weights are recalibrated. Each generated
    ## sample is tested for falling into target +/- tolerance for all
    ## features, in which case it is yielded.
    ##
    ## self.tweak_factor controls the speed of recalibration of weights
    def targeted_samples(self):
        means=None

        counter = 0
        while True:
            fstats = FeatureStatistics()
            for i in range(self.samples_per_round):
                counter+=1
                sample = self.sample()

                ## evaluate all features
                values = { k:self.features[k].eval(sample) for k in self.features }

                ## record the features
                fstats.record_features( self.features, values )

                if self.is_good_sample( self.features, values ):
                    yield sample

            last_means=means
            means = fstats.means()

            # modify weight of each targeted feature
            for fid,f in self.features.items():
                new_weight = 1
                try:
                    new_weight = f.weight + self.tweak_factor*( f.target - means[fid] )
                    self._model.set_feature_weight( new_weight, fid )
                except:
                    pass

            self.requires_reinitialization()

    def targeted_sample( self ):
        try:
            assert( self._targeted_samples != None )
        except:
            self._targeted_samples = self.targeted_samples()

        return next(self._targeted_samples)
