#!/usr/bin/env python3

# -----------------------------
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

## @file
#
#  @brief The infrared high-level Python interface

## @package infrared.infrared
#  @copydoc infrared.py

import os
import random
import re
import sys
import inspect
import math
import subprocess
import copy

from . import libinfrared

from treedecomp import TreeDecompositionFactory
from treedecomp import seed as tdseed

from abc import ABC, abstractmethod


def seed(seed = None):
    """Seed random number generator of libinfrared and treedecomp

    This seeds the RNG of lib infrared (C++ side) and as well
    the random number generator used by randomization in the TreeDecomposition,
    both with the same number.
    It does NOT seed Python's global RNG.

    Args:
        seed: integer used as seed

    Note:
        Without argument or seed==None, use Pythons built-in random.seed() to
        generate a seed.
    """

    if seed is None:
        random.seed()
        seed = random.randint(0, 2**31)

    libinfrared.seed(seed)
    tdseed(seed)


#  exception to signal inconsistency, when consistency would be required
class ConsistencyError(RuntimeError):
    def __init__(self, arg):
        """Construct with description

        Args:
            arg: description of error
        """

        self.args = [arg]


###########
# classes to support different algebras;
# switch between optimization and sampling
#

class EvaluationAlgebra(ABC):
    """Algebra for evaluating a constraint network
    """

    @staticmethod
    @abstractmethod
    def function(weighted_function):
        """Translate weighted function to suitable libinfrared function
            Args:
                weighted_function: an object of class WeightedFunction
        """
        pass

    def interpret(self, value):
        return value

class PFFunctionAdapter(libinfrared.Function):
    """Adapt function for partition function algebra

    @see PFEvaluationAlgebra
    """

    def __init__(self, wf):
        super().__init__(wf.vars())
        self._wf = wf

    def __call__(self, a):
        return math.exp(self._wf.weight * self._wf.value(a))

    def __str__(self):
        return f"<PFFunctionAdapter of {self._wf}>"

class PFEvaluationAlgebra(EvaluationAlgebra):
    """Partition function algebra for sampling

    Controls the translation of WeightedFunction to libinfrared.Function
    in the case of partition function computation / sampling
    """

    def function(self, weighted_function):
        return PFFunctionAdapter(weighted_function)

class ArcticFunctionAdapter(libinfrared.IntFunction):
    """Adapt function for maximization algebra

    @see ArcticEvaluationAlgebra
    """

    def __init__(self, wf, scale):
        super().__init__(wf.vars())
        self._wf = wf
        self._scale = scale

    def __call__(self, a):
        return int(self._wf.weight * self._wf.value(a) * self._scale)

    def __str__(self):
        return f"<ArcticFunctionAdapter of {self._wf}>"

class ArcticEvaluationAlgebra(EvaluationAlgebra):
    """Maximization algebra for optimization

    Controls the translation of WeightedFunction to libinfrared.Function
    in the case of optimization
    """

    def __init__(self, scale):
        self._scale = scale

    def function(self, weighted_function):
        return ArcticFunctionAdapter(weighted_function, self._scale)

    def interpret(self, value):
        return value/self._scale

class WeightedFunction:
    """Function of a constraint network

    Has features weight and variables; their value depends
    on assignment to its variables

    Is 'translated' to a Function of infrared by FunctionAdaptors
    """

    def __init__(self, variables):
        self._vars = variables
        self._weight = 0

    @abstractmethod
    def value(self):
        pass

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, weight):
        self._weight = weight

    def vars(self):
        return self._vars

def _generic_def_function_class(classname, init, value, module="__main__",
                                parentclass=WeightedFunction,
                                valuefunname="value"):
    """Create a class of infrared weighted functions or constraints"""
    def _init(self, *args, **kwargs):
        if "__direct_super__" in kwargs:
            del kwargs["__direct_super__"]
            super(self.__class__, self).__init__(*args, **kwargs)
            return

        variables = init(*args, **kwargs)
        super(self.__class__, self).__init__(variables)

        siginit = inspect.signature(init)
        sigvalue = inspect.signature(value)
        for i, kw in zip(range(len(args)), siginit.parameters):
            kwargs[kw] = args[i]
        self._args = {k: kwargs[k] for k in kwargs if k in sigvalue.parameters}

    def _value(self, a):
        a = a.values()
        params = [a[var] for var in self.vars()]
        return value(*params, **self._args)

    def _str(self):
        return '{} on {}'.format(self.__class__, self.vars())

    def _copy(self):
        cp = (type(self))(self.vars(), __direct_super__=True)
        cp.__dict__.update(self.__dict__)
        return cp

    def _deepcopy(self, memo):
        cp = (type(self))(self.vars(), __direct_super__=True)
        cp.__dict__.update(self.__dict__)
        return cp

    newclass = type(classname, (parentclass,),
                    {
        "__init__": _init,
        valuefunname: _value,
        "__str__": _str,
        "__copy__": _copy,
        "__deepcopy__": _deepcopy,
        "name": lambda self: classname
    })

    sys.modules[module].__dict__[classname] = newclass

def def_function_class(classname, init, value, module="__main__"):
    """
    Define a function class (of type WeightedFunction)

    Defines a new class with name classname (by default, in the main namespace)

    Objects of the defined class can be used in the construction of the
    infrared constraint model. Note that `init` defines the dependencies in
    the order of the (positional) arguments of `value`.

    The value function can depend on arguments to the function init,
    which will be automatically stored in the class and passed on.

    Args:
        classname: name of the class to be defined
        init:      init function of generated class
        value:     value function of generated class
        module:    module where class is generated

    Example:

    @code{.py}
    def_function_class('GCCont', lambda i: [i], lambda x: 2 if x == 1 or x == 2 else 0)
    @endcode

    More complex examples can be found in the accompanying Jupyter tutorials.

    """
    _generic_def_function_class(
        classname, init, value, module, WeightedFunction, "value")


def def_constraint_class(classname, init, value, module="__main__"):
    """
    Define a Constraint class

    Defines a new class with name `classname` (by default, in the main namespace)

    Objects of the defined class can be used in the construction of the
    infrared constraint model. Note that `init` defines the dependencies in
    the order of the (positional) arguments of `value`.

    The value function can depend on arguments to the function init,
    which will be automatically stored in the class and passed on.

    Args:
        classname: name of the class to be defined
        init:      init function of generated class
        value:     value function of generated class
        module:    module where class is generated

    Example:
    @code{.py}
    _bpcomp_tab = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
    def_constraint_class('BPComp', lambda i, j: [i, j],
        lambda x, y: (x, y) in _bpcomp_tab)
    @endcode

    More complex examples can be found in the accompanying Jupyter tutorials.
    """
    _generic_def_function_class(
        classname, init, value, module, libinfrared.Constraint, "__call__")

# -----
# constraint: restrict domain to specific values
class ValueIn(libinfrared.Constraint):
    """
    Constrain variable to have a value from a specified set

    ```
    ValueIn(i, [0,2,3])
    ```
    """
def_constraint_class('ValueIn', lambda i, values: [i],
                     lambda x,values: x in values,
                     module=__name__)

# support special functionality of propagation to domain and entailment
# check when adding this constraint
def _domain_constraint_on_add(self, model):
    i = self.vars()[0]
    values = self._args['values']
    model.restrict_domains(i,(min(values),max(values)))
ValueIn.on_add = _domain_constraint_on_add
def _domain_constraint_entailed(self, model):
    i = self.vars()[0]
    values = self._args['values']
    domain = model.domains[i]
    is_entailed = all( x in values for x in range(domain.lb(), domain.ub()+1) )
    return is_entailed
ValueIn.entailed = _domain_constraint_entailed


class Model:
    """A constraint model

    Describes a constraint model.

    Allows to compose a constraint model by adding variables,
    constraints, and functions and restricting domains.

    Automatically generated features for function groups.

    Provides further information on the model.
    """

    def __init__(self, number=None, domain=None, name='X'):
        """Init model.

        If number is not None, calls add_variables with the
        given parameters after initialization.

        Args:
            number: if not None, number of variables
            domain: domains of variables
            name: name of variable series. Defaults to `X`
        """
        self._constraints = []
        self._functions = dict()
        self._domains = dict()

        self._features = dict()

        if number is not None:
            self.add_variables(number, domain, name)

    def add_variables(self, number, domain, name='X'):
        """Add variable domains

        Args:
            number: number of variables
            domain: domain size; defines the domain values as 0..domain-1
            name: assign a name to the variable(s)
        """

        if name not in self._domains:
            self._domains[name] = []

        self._domains[name].extend(
            [libinfrared.FiniteDomain(domain) for _ in range(number)])

    def restrict_domains(self, vars, domain):
        """Restrict the domain of a variable

        Args:
            vars: variable or list of variables, each specified by
        (name,index); or simply index, then addressing ('X',index)
            domain: the domain

        Note:
            the domain bounds are intersected with the original domain
        """
        newdom = libinfrared.FiniteDomain(domain)
        if type(vars) != list:
            vars = [vars]
        for v in vars:
            name, i = v if type(v) == tuple else ('X', v)

            newlb = max(self._domains[name][i].lb(), newdom.lb())
            newub = min(self._domains[name][i].ub(), newdom.ub())

            self._domains[name][i] = libinfrared.FiniteDomain(newlb,newub)

    def add_constraints(self, constraints):
        """Add constraints to the model
        Args:
            constraints: an iterable of constraints or a single constraint

        Note:
            supports optimizations via on_add and entailed methods of
        added constraints; see ValueIn
        """
        if hasattr(constraints, '__iter__'):
            constraints = list(constraints)
        else:
            constraints = [constraints]

        for constraint in constraints:
            if hasattr(constraint,"on_add"):
                constraint.on_add(self)

            if hasattr(constraint,"entailed"):
                if constraint.entailed(self):
                    continue

            self._constraints.append(constraint)

    def add_functions(self, functions, group='base'):
        """Add functions to the model
        Args:
            functions: [const] an iterable of constraints or a single constraint
            group: indentifier of function group. Defaults to `base`

        Note:
            deep copies the input functions
        """

        if group not in self._functions:
            self._functions[group] = []

        # reserve auto feature entry for group or invalidate cached feature
        self._features[group] = None

        if hasattr(functions, '__iter__'):
            functions = list(functions)
        else:
            functions = [functions]

        self._functions[group].extend(copy.deepcopy(functions))

    def num_named_variables(self, name):
        """Number of variables
        Args:
            name: name of the variables series
        Returns:
            number of variables of the given series
        """
        if name not in self._domains:
            return 0
        return len(self._domains[name])

    def has_empty_domains(self):
        """Check inconsistency due to empty domains
        Returns:
            whether model has empty domains
        """
        return any(dom.empty() for dom in self.domains)

    @property
    def num_variables(self):
        """Number of all variables
        Returns:
            number of all variables
        """
        return sum(len(self._domains[name]) for name in self._domains)

    @property
    def constraints(self):
        """Constraints of the model
        Returns:
            specification of the model's functions
        """
        return self._constraints

    @property
    def functions(self):
        """All functions of the model
        Returns:
            list of all functions
        """
        fs = []
        for k in self._functions:
            fs.extend(self._functions[k])
        return fs

    @property
    def domains(self):
        """Domains of the model
        Returns:
            list of all domain descriptions
        """
        doms = []
        for k in sorted(self._domains.keys()):
            doms.extend(self._domains[k])
        return doms

    def _automatic_feature(self, group):
        """Automatic feature for function group.
        Feature with value that is derived as sum of the feature functions

        Returns:
            the feature
        """
        def eval_fun(sample):
            return sum(f.value(sample) for f in self._functions[group])

        return Feature(group, eval_fun, group=group)

    def add_feature(self, name, group, eval_fun):
        """Add a (custom) feature
        Args:
            name: name of the feature
            group: one or several groups of feature controlled functions
            eval_fun: function to evaluate the feature at an assignment
        """
        self._features[name] = Feature(name, eval_fun, group=group)

    @property
    def features(self):
        """Features of the model
        Returns:
            dictionary of features
        """
        for name in self._features:
            if self._features[name] is None:
                self._features[name] = self._automatic_feature(name)

        return self._features

    def eval_feature(self, assignment, name):
        """Evaluate named feature at assignment
        Args:
            assignment: the assignment
            name: name of the feature
        Returns:
            value of the feature
        """
        return self.features[name].eval(assignment)

    def set_feature_weight(self, weight, name):
        """Set the weight of a feature and its function group.

        This method sets the weight for the Feature itself and in
        all the functions (currently!) in its group.

        The method should be called after all functions of its group are
        defined. Otherwise, weights of the feature and its functions get out
        of sync (but can be re-synced by another call to this method.)

        Args:
            weight: the weight
            name: the feature name
        """
        self.features[name]._weight = weight
        groups = self.features[name].group
        if type(groups) != list:
            groups = [groups]
        for group in groups:
            for f in self._functions[group]:
                f.weight = weight

    def idx(self, variables):
        """Raw indices of named variables
        Args:
            variables: single variable or list of variables; variables can be named; default
        name 'X'
        Returns:
            list of (internal) variable indices
        """
        def convert(var):
            try:
                (name, idx) = var
            except TypeError:
                return var
            offset = 0
            for k in sorted(self._domains.keys()):
                if k == name:
                    break
                offset += len(self._domains[k])
            return offset + idx

        if type(variables) != list:
            variables = [variables]

        variables = [convert(var) for var in variables]
        return variables

    def dependencies(self, non_redundant=True):
        """Dependencies due to constraints and functions
        Args:
            non_redundant: whether the dependency list is made non-redundant. Defaults to `True`
        Returns:
            list of lists of indices of variables that depend on each
        other either through functions or constraints
        """
        deps = [x.vars() for x in self.functions + self.constraints]
        if non_redundant:
            deps = self._remove_subsumed_dependencies(deps)
            deps = list(set(map(tuple, deps)))

        return deps

    def bindependencies(self, non_redundant=True):
        return self._expand_to_cliques(self.dependencies(non_redundant))

    @staticmethod
    def _expand_to_cliques(dependencies):
        """Expand non-binary dependencies to cliques of binary deps
        Args:
            dependencies: list of dependencies
        Returns:
            list of binary dependencies
        """
        import itertools
        bindeps = list()
        for d in dependencies:
            bindeps.extend(itertools.combinations(d, 2))
        return bindeps

    @staticmethod
    def _remove_subsumed_dependencies(deps):
        """Removes redundant dependencies that are subsumed by others
        Removes all dependencies that are subsumed by other dependencies
        in deps

        Args:
            deps: list of dependencies (where a dependency is a list of
                                          indices)

        Returns:
            pruned list of dependencies
        """
        def sublist(xs, ys):
            return all(x in ys for x in xs)

        def subsumed(dep, deps):
            return any(len(dep) < len(dep2) and sublist(dep, dep2)
                       for dep2 in deps)
        deps = [dep for dep in deps if not subsumed(dep, deps)]
        return deps

    def write_graph(self, out, non_redundant=True):
        """Write dependency graph.
        Writes the dependency graph to file in dot format;
        hyper-edges are expanded to cliques.

        Args:
            out: a filehandle of name of the target file
        """
        if type(out) == str:
            out = open(out, 'w')

        out.write("graph G {\n\n")

        for name in sorted(self._domains.keys()):
            for idx, domain in enumerate(self._domains[name]):
                label = name+str(idx)
                out.write("\tvar{} [label=\"{}\"];\n".format(idx, label))
        out.write("\n\n")

        for dep in self._expand_to_cliques(self.dependencies(non_redundant)):
            edgelabel = ''
            x, y = dep
            out.write(
                "\tvar{} -- var{}  [label=\"{}\"];\n".format(x, y,
                                                             edgelabel))

        out.write("\n}\n")

    def connected_components(self):
        """Connected components of the model's dependency graph
        Returns:
            a list of sets of the connected components
        """
        numnodes, edges = self.num_variables, self.bindependencies()
        def adjacency_list():
            al = [[] for _ in range(numnodes)]
            for d in edges:
                for x in d:
                    al[x].extend([y for y in d if x!=y])
            return [sorted(set(xs)) for xs in al]

        al = adjacency_list()
        marked = [False] * numnodes
        def component(x):
            if marked[x]:
                return []
            marked[x] = True
            c = [x]
            for y in al[x]:
                c.extend(component(y))
            return c
        return [set(component(x)) for x in range(numnodes)
                if not marked[x]]

class ClusterTreeBase:
    """Cluster tree base class

    This class provides functionality to construct and populate C++
    cluster trees with constraints and functions.

    It is used as mix-in class for the specialized cluster tree classes,
    which wrap interface the C++/libinfrared cluster tree classes.
    """

    def __init__(self, model, td, EvaluationAlgebra):
        self._model = model
        self._td = td
        self._EA = EvaluationAlgebra

        self._bagsets = list(map(set, td.get_bags()))

        self.construct_cluster_tree(model.domains, td)

    def evaluate(self):
        return self._EA.interpret(self._ct.evaluate())

    def is_consistent(self):
        """
        Returns:
            whether the cluster tree is consistent
        """
        return self._ct.is_consistent()

    def construct_cluster_tree(self, domains, td):
        """Construct the cluster tree object of the C++ engine

        Args:
            domains: description of the domains
            td: tree decomposition
        """
        bagconstraints, bagfunctions = self.get_bag_assignments()

        # keep record of all non-root nodes
        children = set()

        for bagidx in td.toposorted_bag_indices():
            if bagidx not in children:  # --> bagidx is a root
                # perform DFS of subtree of bagidx:
                #   add cluster to cluster tree
                #    and populate them with functions and constraints
                # stack entries: parent cluster index, bag index
                stack = [(None, bagidx)]
                while stack:
                    (p, i) = stack.pop()
                    bagvars = sorted(list(self._bagsets[i]))

                    if p is None:
                        cluster = self._ct.add_root_cluster(bagvars)
                    else:
                        cluster = self._ct.add_child_cluster(p, bagvars)

                    for x in bagconstraints[i]:
                        self._ct.add_constraint(cluster, x)

                    for x in bagfunctions[i]:
                        self._ct.add_function(cluster, self._EA.function(x))

                    for j in td.adj[i]:
                        children.add(j)
                        stack.append((cluster, j))

    def get_bag_assignments(self):
        """Get assignments of functions and constraints to the bags

        straightforward non-redundant assignment of constraints and functions,
        each to some bag that contains all of their variables

        assumes constraints and functions specified in self._model
        """
        bagconstraints = self.assign_to_bags(self._model.constraints)
        bagfunctions = self.assign_to_bags(self._model.functions)
        return (bagconstraints, bagfunctions)


    def find_all_bags(self, bvars):
        """Get the indices of all bags that contain a set of variables
        Args:
            bvars: the set of variables
        Returns:
            list of indices of the bags that contain bvars
        """
        return [i for i, bag in enumerate(self._bagsets)
                if all(x in bag for x in bvars)]

    def find_bag(self, bvars):
        """Find a bag that contains a set of variables
        Args:
            bvars: the set of variables
        Returns:
            index of first bag that contains bvars (or None if there is none)
        """
        bags = self.find_all_bags(bvars)
        if len(bags) > 0:
            return bags[0]
        else:
            return None

    def assign_to_bags(self, constraints):
        """assign constraints or functions to bags

        Assigns such that each constraint/function is assigned to exactly one bag
        that contains its dependencies
        Args:
            constraints: list of constraints
        Returns:
            list where constraints are placed at corresponding bag indices

        @pre for each constraint/function there is one bag that contains its
            dependencies; otherwise the constraint/function is not assigned
        """
        bagconstraints = {i: [] for i in range(len(self._bagsets))}
        for cons in constraints:
            bagconstraints[self.find_bag(cons.vars())].append(cons)
        return bagconstraints

    def assign_to_all_bags(self, constraints):
        """ assign constraints or functions to all possible bags

        Assigns constraint/function is to each bags that contains its
        dependencies.
        Args:
            constraints: list of constraints
        Returns:
            list where constraints are placed at corresponding bag indices
        """
        bagconstraints = {i: [] for i in range(len(self._bagsets))}
        for cons in constraints:
            for bidx in self.find_all_bags(cons.vars()):
                bagconstraints[bidx].append(cons)
        return bagconstraints

class ArcticClusterTree(ClusterTreeBase):
    """
    Cluster tree for maximization (in arctic semiring)

    @see infrared.infrared.ArcticOptimizer
    """
    def __init__(self, model, td, scale = 100):
        self._ct = libinfrared.ArcticClusterTree(model.domains)
        super().__init__(model, td, ArcticEvaluationAlgebra(scale))

    def optimize(self):
        return self._ct.optimize()


class PFClusterTree(ClusterTreeBase):
    """
    Cluster tree for partition function calculation and sampling
    """
    def __init__(self, model, td):
        self._ct = libinfrared.PFClusterTree(model.domains)
        super().__init__(model, td, PFEvaluationAlgebra())

    def sample(self):
        return self._ct.sample()


class Feature:
    """Feature in multi-dimensional Boltzmann sampling

    A feature defines a (partial) evaluation/score of a sample, which
    can be targeted due to a dedicated weight. It defines a method
    value, which determines the feature's value for a sample.

    Moreover, a feature controls one or several groups of functions.

    Features should belong to exactly one Model. Features have weights,
    but their weights have to be synchronized with their functions.
    The weight must be changed only by the model in sync
    with the corresponding functions (set_feature_weight). For this reason,
    weight is defined as a read-only property.
    """

    def __init__(self, identifier, eval_function, group=[]):
        """Construct feature

        Args:
            ifentifier:
            eval_funciton:
            group: defaults to empty list []
        """
        self._identifier = identifier
        self._eval_function = eval_function
        self._weight = 0
        self._group = group

    @property
    def identifier(self):
        return self._identifier

    def eval(self, sample):
        """Evaluate feature for given sample and weight"""
        return self._eval_function(sample)

    @property
    def group(self):
        return self._group

    @property
    def weight(self):
        return self._weight


class FeatureStatistics:
    """Keeping statistics on features

    This class allows recording values of multiple features for a
    series of samples; it can be queried for mean and standard
    deviation (or the entire distribution) of each recorded feature.
    """

    def __init__(self, keep=False):
        """Construct

        Args:
            keep: Keep all recorded features in memory. Defaults to False
        """
        self.keep = keep

        self.identifier = dict()
        self.count = dict()
        self.sums = dict()
        self.sqsums = dict()
        self.features = dict()

    def record_features(self, features, values):
        """Record feature values

        Args:
            features: a dictionary of features
            values: a corresponding dictionary of the feature values
        Returns:
            pair of feature id string and value
        """
        if not type(features) is dict:
            features = {'dummy': features}
            values = {'dummy': values}

        for k in features:
            value = values[k]
            fid = features[k].identifier

            if fid not in self.count:
                self.identifier[fid] = features[k].identifier
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

    def empty(self):
        """Check whether any features have been recorded

        Returns:
            whether empty (since no feature has been recorded)
        """
        return not self.count

    def means(self):
        """Means of recorded features

        Returns:
            dictionary of means  (at feature ids as keys)
        """
        return {k: self.sums[k]/self.count[k] for k in self.sums}

    def variances(self):
        """Variances of recorded features
        Returns:
            dictionary of variances (at feature ids as keys)
        """
        means = self.means()
        return {k: self.sqsums[k]/self.count[k] - means[k]**2
                for k in self.sqsums}

    def stds(self):
        """Standard deviations of recorded features
        Returns:
            dictionary of standard deviations (at feature ids as keys)
        """
        return {k: val**0.5 for k, val in self.variances().items()}

    def report(self):
        """Report features to standard output
        """
        means = self.means()
        stds = self.stds()
        return ' '.join(["{}={:3.2f} +/-{:3.2f}".format(self.identifier[fid],
                                                        means[fid],
                                                        stds[fid])
                         for fid in self.count])


class EngineBase(ABC):
    """Abstract base class for samplers and optimizers
    """

    #  Construct from model
    #  model the constraint model
    def __init__(self, model, td_factory=TreeDecompositionFactory(),
                 lazy=True):
        """Construct with model and optional td_factory

        Args:
            model: [const] Constraint network model
            td_factory: Factory for tree decomposition
            lazy: delay construction of the data structures until required

        Note:
            the model is deepcopied such that it won't be modified and/or
        later modifications of the model don't have effect on the sampler
        """

        self._model = copy.deepcopy(model)
        self._td_factory = td_factory

        self.requires_reinitialization()
        if not lazy:
            self.setup_engine()

    #  flag that engine requires setup
    def requires_reinitialization(self):
        self._td = None
        self._ct = None

    @property
    def model(self):
        return self._model

    @property
    def td(self):
        self.setup_engine(skip_ct=True)
        return self._td

    @property
    def ct(self):
        return self._ct

    def setup_engine(self, *, skip_ct=False):
        """Sets up the constraint model / cluster tree sampling
            engine

        Note:
            do nothing, if the engine was already initialized before
        Args:
            skip_ct: skip the potentially expensive construction and
        evaluation of the cluster tree
        """
        # immediately return if ct exists and is not None
        if self._ct is not None:
            return

        if self._td is None:
            self._td = self._td_factory.create(
                self._model.num_variables, self._model.dependencies())

        if not skip_ct:
            if self._model.has_empty_domains():
                raise ConsistencyError("Model has empty domains")

            self._ct = self.gen_cluster_tree()

            if not self._ct.is_consistent():
                raise ConsistencyError("Inconsistent constraint model")

    def evaluate(self):
        """Evaluates the cluster tree
        Returns:
            partition function
        Note:
            Evaluation is a potentially (depending on the treewidth)
        costly operation.
        The method does not re-evaluate the tree if this was already done
        """
        self.setup_engine()
        return self._ct.evaluate()

    def is_consistent(self):
        try:
            self.setup_engine()
        except ConsistencyError:
            return False
        return self._ct.is_consistent()

    def plot_td(self, filename, to="pdf"):
        """Plot the tree decomposition to pdf file
        filename write to filename
        to target format, support conversion to "pdf" or "png".
        Anything else writes graphviz dot format
        """
        conversions = {
            "pdf": dotfile_to_pdf,
            "png": dotfile_to_png
        }

        if to in conversions:
            filename = re.sub(f".{to}$", "", filename)+".dot"

        self.setup_engine(skip_ct=True)
        with open(filename, "w") as dot:
            self._td.writeTD(dot)

        if to in conversions:
            conversions[to](filename)
            os.remove(filename)

    def treewidth(self):
        """Get tree width
        Returns:
            tree width
        """
        self.setup_engine(skip_ct=True)
        return self._td.treewidth()


    @abstractmethod
    def gen_cluster_tree(self):
        """Generate the populated cluster tree
        td tree decomposition
        Returns: cluster tree
        """
        pass


class ArcticOptimizer(EngineBase):
    """Maximizing optimizer (based on arctic algebra)
    """

    def __init__(self, model, td_factory=TreeDecompositionFactory(),
                 lazy=True):
        """ Construct
        @see EngineBase
        """
        super().__init__(model, td_factory, lazy)

    def optimize(self):
        """Optimal assignment
        Returns:
            one optimal assignment
        """
        self.setup_engine()
        return self._ct.optimize()

    def gen_cluster_tree(self):
        """!  Suitable cluster tree
        Returns:s arctic cluster tree for the model
        """
        return ArcticClusterTree(self._model, td=self._td)


## Alias for default optimizer ArcticOptimizer
Optimizer = ArcticOptimizer


class BoltzmannSampler(EngineBase):
    """Boltzmann sampler
    """

    def __init__(self, model, td_factory=TreeDecompositionFactory(),
                 lazy=True):
        """Construct
        Args:
            model: Infrared Model
            td_factory: Tree decomposition policy. Defaults to TreeDecompositionFactory
            lazy: Defaults to True
        See also:
            EngineBase
        """
        super().__init__(model, td_factory, lazy)

    def sample(self):
        """Generate sample

        If the cluster tree was not evaluated (or consistency checked) before,
        it will be evaluated once on-demand.
        Returns:
            a raw sample

        Note:
            raises exception ConsistencyError if the model is inconsistent.
        """
        self.setup_engine()
        return self._ct._ct.sample()

    def samples(self):
        """Sample generator
        """
        while(True):
            yield self.sample()

    def gen_cluster_tree(self):
        """Suitable cluster tree
        Returns:
            PFClusterTree for the model
        """
        return PFClusterTree(self._model, td=self._td)


class MultiDimensionalBoltzmannSampler(BoltzmannSampler):
    """ Multi-dimensional Boltzmann sampler
    """

    def __init__(self, model, td_factory=TreeDecompositionFactory(),
                 lazy=True):
        """Construct with model and optional td_factory

        Args:
            model: [const] Constraint network Model
            td_factory: Factory for tree decomposition
            lazy: delay construction of the data structures until required. Defaults to True

        See also:
            BoltzmannSampler.__init__()
        """
        super().__init__(model, td_factory, lazy)

        # parameters controlling the mdbs procedure
        self.samples_per_round = 1000
        self.tweak_factor = 0.01

    def is_good_sample(self, features, values):
        """whether the sample is of good quality

        checks whether the sample approximately meets the targets;
        check only the targeted features (which have value, target and tolerance)

        Args:
            features: dictionary of features
        """
        ret = True
        for k, f in features.items():
            try:
                if abs(values[k] - f.target) > f.tolerance:
                    ret = False
                    break
            except AttributeError:
                pass
        return ret

    def set_target(self, target, tolerance, featureid):
        """Set target of a feature

        Args:
            target: the target value
            tolerance: the tolerance (as absolute difference) to the target
            fetureid: id of the feature
        """
        f = self._model.features[featureid]
        f.target = target
        f.tolerance = tolerance

    def targeted_samples(self):
        """Generator of targeted samples

        Performs multi-dimensional Boltzmann sampling: every
        self.samples_per_round many samples, the feature means are
        estimated and the weights are recalibrated. Each generated
        sample is tested for falling into target +/- tolerance for all
        features, in which case it is yielded.

        self.tweak_factor controls the speed of recalibration of weights
        """
        features = self._model.features

        means = None

        counter = 0
        while True:
            fstats = FeatureStatistics()
            for i in range(self.samples_per_round):
                counter += 1
                sample = self.sample()

                # evaluate all features
                values = {k: features[k].eval(sample) for k in features}

                # record the features
                fstats.record_features(features, values)

                if self.is_good_sample(features, values):
                    yield sample

            means = fstats.means()

            # modify weight of each targeted feature
            for fid, f in features.items():
                new_weight = 1
                try:
                    new_weight = (f.weight +
                                  self.tweak_factor * (f.target - means[fid]))
                    self._model.set_feature_weight(new_weight, fid)
                except AttributeError:
                    pass

            self.requires_reinitialization()

    def targeted_sample(self):
        if not hasattr(self,'_targeted_samples'):
            self._targeted_samples = self.targeted_samples()

        return next(self._targeted_samples)


## Alias for MultiDimensionalBoltzmannSampler
Sampler = MultiDimensionalBoltzmannSampler


def dotfile_to_tgt(graphfile, tgt, outfile=None):
    """Convert dot graph file format to png/pdf

    The graph is plotted and written to a pdf or png file by calling
    graphviz's dot tool on the dot file.

    Args:
        graphfile: file of graph in dot format
        tgt: target format ('png' or 'pdf')
        outfile: optional output filename

    Note:
        Unless specified, the output filename is composed from the input file name
        and the target format.
    """

    if outfile is None:
        outfile = re.sub(r".dot$", "", graphfile)+"."+tgt
    subprocess.check_output(["dot", f"-T{tgt}", "-o", outfile, graphfile])


def dotfile_to_pdf(graphfile, outfile=None):
    """Convert dot graph file to pdf

    @see dotfile_to_tgt
    """

    dotfile_to_tgt(graphfile, "pdf", outfile)


def dotfile_to_png(graphfile, outfile=None):
    """Convert dot graph file to png

    @see dotfile_to_tgt
    """

    dotfile_to_tgt(graphfile, "png", outfile)

