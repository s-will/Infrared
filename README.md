# Infrared

[TOC]

Infrared is a generic C++/Python hybrid library for efficient
(fixed-parameter tractable) Boltzmann sampling.

## Disclaimer and license

Infrared is free software. It was part of project [RNARedPrint](https://github.com/yannponty/RNARedPrint), then separated as a stand alone project.
Note that the system is in active development and is likely to still undergo significant
changes. It is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details
[<http://www.gnu.org/licenses/>].

## Main features of Infrared

Infrared efficiently solves a broad class of sampling and optimization
problems; originally targeted to bioinformatics problems, especially
RNA sequence design. Such problems can be specified in a declarative, compositional
style as (evaluated) constraint models using the Python high-level
interface of Infrared.

Encapsulated by the Python interface, the framework provides a fast and
flexible C++ engine that evaluates constraint networks consisting of
variables, multi-ary functions, and multi-ary constraints.  This evaluation
is efficient depending on the complexity, measured as tree-width, of the
network. 

Application specific functions and constraints are easily defined in Python
(or in C++). The evaluation is performed efficiently using cluster
tree elimination following a (hyper-)tree decomposition of the
dependencies (due to functions and constraints). Interpreting the
evaluations as partition functions, the system supports sampling of
variable assignments from the corresponding Boltzmann distribution.
Finally, Infrared implements a generic multi-dimensional Boltzmann
sampling strategy to target specific feature values. Such
functionality is made conveniently available via general Python
classes. In particular, the interface was defined to allow the
straightforward and declarative specification of the constraint model.

### Characteristics in brief

Keep it simple, but flexible on the C++ side

- variables are indexed from 0..n-1; using index type int

- variable domains are consecutive lb...ub; value type is int

- function value type templated

- function values are combined according to policy, which supports evaluation
  algebras

- support sparse, memory-optimized data structures for intermediary
  results to enable larger problem instances

Provide a flexible interface and added functionality on the Python side

- support constraint modeling syntax to declarativly and compositionally
  describe constraint problems

- support definintion of 'custom' constraints and functions in Python

- support efficient evaluation of models for optimization and Boltzmann sampling

## Installation

There are two ways to install Infrared. The first option is using conda.
Infrared is depolyed on conda-forge channel.
Users can skip the first line command if it's already done.

```
conda config --add channels conda-forge 
conda install infrared
```


For users who don't want to use conda, Infrared can also be installed with standard pip install.
The compilation requires [pybind11](https://github.com/pybind/pybind11) (at least v2.4) and CMake.
Note that pybind11 installed from standard ubuntu18.04 APT is outdated and the variable `PYBIND11_GLOBAL_SDIST` should be set for version before v2.6.0 as 

```
PYBIND11_GLOBAL_SDIST=1 python3 -m pip install https://github.com/pybind/pybind11/archive/master.zip
python3 -m pip install .
```

## Usage

We provide a tutorial as introduction to the Python high-level interface in
 the jupyter notebook [Doc/rnadesign-tutorial.ipynb](https://gitlab.inria.fr/amibio/Infrared/-/tree/master/Doc/rnadesign-tutorial.ipynb) and further
 examples notebooks in the [Doc subdirectory](https://gitlab.inria.fr/amibio/Infrared/-/tree/master/Doc).

## API documentation

The Infrared API is documented using doxygen comments, such that
html documentation can be generated (in Doc/html) by doxygen.

A further entry-point to using the library in novel sampling applications
is provided by the code of RNARedPrint 2.

More information on the library architecture is provided below.

## Infrared architecture and background

The system was build to separate the core "Infrared" from applications
like [RNARedprint 2](https://gitlab.inria.fr/amibio/RNARedPrint) or
[RNAPOND](https://gitlab.inria.fr/amibio/RNAPOND).


### infrared / libinfrared --- Python high level interface

Using pylib11, we expose classes of the C++ library to
Python, such that cluster trees can be constructed and populated with
exposed C++ constraints/functions and/or derived Python
constraints/functions. 
Python programs  using the infrared library,
typically import module infrared, create an instance of Model and
populate it with variables (specifying their finite domains), constraints
and functions. The model automatically generates features from the
functions; in order to generate and control several features, the functions
can be assigned to function groups. Moreover, the user can define
additional features.

The model is passed to Samplers that can generate samples
from a specific Boltzmann distrubition as well as samples targeted at
specific feature values. The latter performed by multi-dimensional
Boltzmann sampling. The module provides access to the Infrared core and
exports the additional base classes 


### ired --- Infrared C++ core

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

By design, the Infrared core is strictly low-level and
domain-agnostic, it only knows finite domain variables, constraints,
and functions; as well as to evaluate partition functions and sample
based on cluster trees for constraint networks. This will allow using
the same core for evaluation (even not necessarily of partition
functions, e.g. instead optimizing energy/fitness) and sampling in
very different domains just by writing domain-specific Python code.
