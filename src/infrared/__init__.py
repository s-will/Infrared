#
# InfraRed -  A generic engine for Boltzmann sampling over constraint networks
# (C) Sebastian Will, 2018
#
# This file is part of the InfraRed source code.
#
# InfraRed provides a generic framework for tree decomposition-based
# Boltzmann sampling over constraint networks
#

## @file __init__.py
# Package infrared

""" @package infrared
Package infrared

Package defining Infrared's Python interface.

Note that many names from submodules infrared and some from libinfrared are raised \
to the package-level namespace.
For example, `infrared.infrared.Model` has the alias `infrared.Model` and `infrared.libinfrared.Constraint` is available as `infrared.Constraint`. These aliases are made explicit as variables.
"""

from . import libinfrared
from . import infrared

# ---- aliases from libinfrared

## @copydoc libinfrared.Function
Function = libinfrared.Function

## @copydoc infrared.libinfrared.Constraint
Constraint = libinfrared.Constraint

## @copydoc infrared.libinfrared.FiniteDomain
FiniteDomain = libinfrared.FiniteDomain

## @copydoc infrared.libinfrared.Assignment
Assignment = libinfrared.Assignment

# ---- aliases from infrared

## @copydoc infrared.infrared.seed
seed = infrared.seed

## @copydoc infrared.infrared.Model
Model = infrared.Model

## @copydoc infrared.infrared.Feature
Feature = infrared.Feature

## @copydoc infrared.infrared.FeatureStatistics
FeatureStatistics = infrared.FeatureStatistics

## @copydoc infrared.infrared.ArcticClusterTree
ArcticClusterTree = infrared.ArcticClusterTree

## @copydoc infrared.infrared.PFClusterTree
PFClusterTree = infrared.PFClusterTree

## @copydoc infrared.infrared.ArcticOptimizer
ArcticOptimizer = infrared.ArcticOptimizer

## @copydoc infrared.infrared.Optimizer
Optimizer = infrared.Optimizer

## @copydoc infrared.infrared.BoltzmannSampler
BoltzmannSampler = infrared.BoltzmannSampler

## @copydoc infrared.infrared.MultiDimensionalBoltzmannSampler
MultiDimensionalBoltzmannSampler = infrared.MultiDimensionalBoltzmannSampler

## @copydoc infrared.infrared.Sampler
Sampler = infrared.Sampler

## @copydoc infrared.infrared.ConsistencyError
ConsistencyError = infrared.ConsistencyError

## @copydoc infrared.infrared.def_function_class
def_function_class = infrared.def_function_class

## @copydoc infrared.infrared.def_constraint_class
def_constraint_class = infrared.def_constraint_class

## @copydoc infrared.infrared.dotfile_to_pdf
dotfile_to_pdf = infrared.dotfile_to_pdf

## @copydoc infrared.infrared.dotfile_to_png
dotfile_to_png = infrared.dotfile_to_png

## @copydoc infrared.infrared.ValueIn
ValueIn = infrared.ValueIn

__version__ = '1.0'
