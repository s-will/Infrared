from .libinfrared import Function, Constraint
from . import libinfrared as libir
from .infrared import seed, ConstraintNetwork, ClusterTree, Feature,\
    FeatureStatistics, BoltzmannSampler, MultiDimensionalBoltzmannSampler,\
    ConsistencyError, def_function_class, def_constraint_class
from . import rna
from .automaton import words_to_accept, words_to_avoid


__version__ = '0.5.1'
