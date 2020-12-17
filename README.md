# Infrared

Infrared is a generic C++/Python hybrid library for efficient
(fixed-parameter tractable) Boltzmann sampling.

## Disclaimer and license

Infrared is free software. It was part of project [RNARedPrint](https://github.com/yannponty/RNARedPrint), then separated as a stand alone project.
Note that the system is still in an early
stage of development and is likely to still undergo significant
changes. It is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details
[<http://www.gnu.org/licenses/>].

## Main features of Infrared

Infrared provides a fast and flexible C++ engine that evaluates a
constraint network consisting of variables, multi-ary functions, and
multi-ary constraints. Functions and constraints are C++ or Python
objects, where new functions and constraints are easily added in C++
or in Python! The evaluation is performed efficiently using cluster
tree elimination following a (hyper-)tree decomposition of the
dependencies (due to functions and constraints).  Interpreting the
evaluations as partition functions, the system supports sampling of
variable assignments from the corresponding Boltzmann distribution.
Finally, Infrared implements a generic multi-dimensional Boltzmann
sampling strategy to target specific feature values. Such
functionality is made conveniently available via general Python
classes.

## Installation

There are two ways to install Infrared. The first option is using conda.
Infrared is depolyed on conda-forge channel.
Users can skip the first line command if it's already done.

```
conda config --add channels conda-forge 
conda install infrared
```


For users who don't want to use conda, Infrared can also be installed with standard pip install.
The compilation requires [boost.graph](https://www.boost.org/), [pybind11](https://github.com/pybind/pybind11) (at least v2.4) and CMake.
Note that pybind11 installed from standard ubuntu18.04 APT is outdated and the variable `PYBIND11_GLOBAL_SDIST` should be set for version before v2.6.0 as 

```
PYBIND11_GLOBAL_SDIST=1 python3 -m pip install https://github.com/pybind/pybind11/archive/master.zip
python3 -m pip install .
```

## Usage

We provide a tutorial as introduction to the Python high-level interface in
 the jupyter notebook ```infrared-rnadesign-tutorial.ipynb```.

## API documentation

The Infrared API is documented using doxygen comments, such that
html documentation can be generated (in Doc/html) by doxygen
```

A further entry-point to using the library in novel sampling applications
is provided by the code of RNARedPrint 2.

More information on the library architecture is provided below.

## Infrared architecture and background

The system is build to separate the core "Infrared" from applications
like Redprint.

### ired --- Infrared core

The core preforms precomputation, i.e. calculation of partition
functions, and Boltzmann sampling based on a given, populated cluster
tree. A /cluster tree/ of a constraint network (of variables,
functions and constraints) corresponds to a (hyper-)tree decomposition
of the dependency graph (induced by the functions/constraints on the
variables). It consists of bags (aka clusters) of variables,
functions, and constraints. The bags are connected by edges to form a
tree (or forest), such that the following properties hold

1) For each variable in the CN, there is one bag that contains the variable.

2) For each variable, the bags that contain this variable form a subtree.

3) For each function, there is exactly one bag that contains
   the function and its variables.

4) For each constraint, there is at least one bag that contains the
   function and its variables. Constraints are only assigned to bags
   that contain all of their variables.

The infrared core (namespace ired) defines the following classes
```
ired::Assignment           an assignment of variables to values
ired::AssignmentIterator   provides iteration over sub-assignments; the work horse of ired
ired::Function<ValueType>  defines functions over variables; abstract
ired::Constraint           defines constraints over variables; abstract
ired::ConstraintNetwork    implements the constraint network (holds vars, funs, and constrs)
ired::Cluster              bag of variables, functions, constraints
ired::ClusterTree          tree of bags, provides evaluation and sampling
```

By design, the Infrared core is strictly low-level and
domain-agnostic, it only knows finite domain variables, constraints,
and functions; as well as to evaluate partition functions and sample
based on cluster trees for constraint networks. This will allow using
the same core for evaluation (even not necessarily of partition
functions, e.g. instead optimizing energy/fitness) and sampling in
very different domains just by writing domain-specific Python code.

### ired::rnadesign --- RNA design specific extension of Infrared

In namespace ired::rnadesign, the system provides domain-specific
constraints and functions for rna design, as such it specializes the
classes ired::Function<double> and ired::Constraint, e.g.
```
ired::rnadesign::ComplConstraint  complementarity constraint
ired::rnadesign::GCControl        function for control of GC-content
ired::rnadesign::BPEnergy         evaluate base pair energy (base pair model)
ired::rnadesign::StackEnergy      evaluate stacking energy (stacking model)
```

### infrared / libinfrared --- Python modules exposing infrared core and extensions

Using boost::python, we expose classes of ired and ired::rnadesign to
Python, such that cluster trees can be constructed and populated with
exposed C++ constraints/functions and/or derived Python
constraints/functions. Python programs like `redprint.py` that want to
use the infrared library, typically import only module infrared, which
provides base classes that can be specialized to generate
application-specific constraint networks and cluster trees for
Infrared. Moreover the module provides classes to generate samples
from a specific Boltzmann distrubition as well as samples targeted at
specific feature values. The latter performed by multi-dimensional
Boltzmann sampling. The module provides access to the Infrared core and
exports the additional base classes 
```
infrared.Feature            derived classes represent features like GC-content, or Turner energy
infrared.FeatureStatistics  gathers statistsics of the feature values for a series of samples
infrared.ConstraintNetwork  derived to construct the constraint network (variables,
                            constraints, functions) of a problem instance 
infrared.TreeDecomposition  derived to construct and represent the tree decomposition
infrared.BoltzmannSampler   wraps the Boltzmann sampling functionality
infrared.MultiDimensionalBoltzmannSampler
                            additionally implements multi-dimensional Boltzmann sampling
```

