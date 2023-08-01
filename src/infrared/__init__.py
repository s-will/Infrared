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

Note that many names from submodules infrared and some from libinfrared are raised
to the package-level namespace.
For example, `infrared.infrared.Model` is available as `infrared.Model`, `infrared.libinfrared.Constraint` as `infrared.Constraint`...
"""

from . import libinfrared
from . import infrared

globals().update(vars(libinfrared))
globals().update({k:v for k,v in vars(infrared).items() if k in infrared._exports})

__version__ = '1.1'
