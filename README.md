# Infrared

Infrared is a generic C++/Python hybrid library for efficient
(fixed-parameter tractable) Boltzmann sampling.

## Disclaimer and license

Infrared is free software. Note that the system is still in an early
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

While the core library is agnostic to the origin of the tree decomposition,
by default it interfaces networkx (optionally, libhtd and TDlib) to perform
the tree decomposition .

## Installation

The software can be compiled and installed (after cloning the
repository) on GNU/Linux or other Unix-like systems as follows.  In
preparation, install libhtd
<https://github.com/mabseher/htd/releases>. Moreover, compilation
requires boost.graph and boost.python. Compile and install Infrared
itself by

```
cd Infrared
./autoreconf -i
./configure --with-htd=$path_to_htd_installation --prefix=$path_to_infrared_installation
make && make install
```

To use tree decomposition based on the Java tool TDlib, please obtain
a copy and set the Java CLASSPATH accordingly.

## RNARedPrint 2

On top of the Infrared library, we re-implemented
RNARedPrint [<https://github.com/yannponty/RNARedPrint>], which is
shipped together with the library and serves as non-trivial
example for the use of the Infrared engine. RNARedprint samples RNA
sequences from a Boltzmann distribution based on the energies of
multiple target structures and GC-content. Our research paper on
RNARedPrint
> Stefan Hammer, Yann Ponty, Wei Wang, Sebastian Will. Fixed-Parameter Tractable Sampling for RNA Design with Multiple Target Structures. Proc. of the 22nd Annual International Conference on Research in Computational Molecular Biology, RECOMB'18. Apr 2018, Paris, France. 2018.

describes many of the ideas that
lead to the development of Infrared and points to the origins in
cluster tree elimination and constraint networks. If you find this
software useful for your own work, please do not forget to cite us.

## Running RNARedPrint 2 from the command line

For running RNARedprint, one furthermore needs an installation of the
ViennaRNA package; make sure that the Python module is found in the
path described by the environment variable PYTHONPATH. PYTHONPATH must
as well point to the libinfrared module from the Infrared
installation, i.e. $path_to_infrared_installation/lib

Redprint reads the multiple RNA target structures from an 'inp'-file, e.g.

```
....((((((..(((((((....))))((((((......))..))))..((((((....((((...)))).....)))))).))).))))))........
..............((((((.....(((.(((((..((((..(((((...((......))....)))))..))))..))...)))))).)))))).....
..((((.((.((.........((((((((.((.(((......))).)).))))))))..)))).)))).(((......................)))...
;
```

A typical call to produce 20 Boltzmann samples with given weights looks like

```
redprint.py test.inp  -n 20 --model bp --gcweight=0.15 --weight 5 --method 0 --turner
```

Moreover, redprint supports targeting specific gc-content and
energies, which is calculated by performing multi-dimensional
Boltzmann sampling (mdbs). Here is a example call
```
time _inst/bin/redprint.py test.inp -n 20 --mdbs --gctarget=70 --gctol 3 \
        --tar -15 --tar -20 --tar -20 --tol 1 --gc --turner -v
```

Further usage information is available by calling `redprint.py --help`
A further tool `redprint_complexity.py` is provided to report tree
widths for two energy models of different complexity (see our research
paper for full details) and plot the dependency graphs and tree
decompositions. This tool provides insight into the run time / space
consumption behavior of redprint on specific instances.

## API documentation

We provide a tutorial as introduction to the Python high-level interface in
 the jupyter notebook ```infrared-rnadesign-tutorial.ipynb```.

The Infrared API is documented using doxygen comments, such that
html documentation can be generated (in Doc/html) by 
```
make doxygen-doc
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

### treedecomp.py --- Generation of tree decompositions

This python module provided by infrared supports the computation of tree
decompositions based on the python framework networkx or (optionally)
external libs (interfaces the Java lib TDlib and the C++-lib libhtd; for
the latter infrared implements a wrapper (module libhtdwrap), which exposes
a specific variant of libthd-tree decomposition to Python.

### redprint.py --- Multi-dimensional Boltzmann sampling for multi-target RNA design

Redprint uses the Infrared library to perform multi-dimensional
Boltzmann sampling of sequences for given multi-target RNA design
instances. For either the base pair or the stacking model, it computes
the dependencies, calls the tree decomposition construction, and
constructs and populates the cluster tree for Infrared. For this
purpose, it specialized base classes of the module infrared. Finally,
it generates samples and optionally reports (statistics over) features
of the sampled sequences, again with the support of the infrared
module.
