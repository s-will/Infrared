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
objects, where new functions and constraints are easily added in C++ or
in Python! The evaluation is performed efficiently using cluster tree
elimination following a (hyper-)tree decomposition of the dependencies
(due to functions and constraints).  Interpreting the evaluations as
partition functions, the system supports sampling of variable
assignments from the corresponding Boltzmann distribution.

While the core library is agnostic to the origin of the tree
decomposition, it comes with interfaces to TDlib and libhtd.

## RNARedPrint 2

On top of the Infrared library, we re-implemented
RNARedPrint [<https://github.com/yannponty/RNARedPrint>], which is
shipped together with the library and serves as first non-trivial
example for the use of the Infrared engine. RNARedprint samples RNA
sequences from a Boltzmann distribution based on the energies of
multiple target structures and GC-content. Our research paper on
RNARedPrint "Stefan Hammer, Yann Ponty, Wei Wang, Sebastian Will.
Fixed-Parameter Tractable Sampling for RNA Design with Multiple Target
Structures. Proc. of RECOMB, 2018." describes many of the ideas that
lead to the development of Infrared and points to the origins in
cluster tree elimination and constraint networks. If you find this
software useful for your own work, please do not forget to cite us.

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

## API documentation

Generate html documentation of the Infrared API by
```
make doxygen-doc
```
Documentation will be generated in Doc/html; find the main page at Doc/html/index.html

## Running Redprint

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

A typical call to produce 20 Boltzmann samples looks like

```
redprint.py test.inp  -n 20 --model bp --gcweight=0.15 -w5 --method 0 --turner
```

Further usage is given by redprint.py --help

```
usage: redprint.py [-h] [--method METHOD] [-n NUMBER] [-v] [--model MODEL]
                   [--turner] [--checkvalid] [--no_redundant_constraints]
                   [--gcweight GCWEIGHT] [-w WEIGHT] [--plot_td]
                   infile

Boltzmann sampling for RNA design with multiple target structures

positional arguments:
  infile                Input file

optional arguments:
  -h, --help            show this help message and exit
  --method METHOD       Method for tree decomposition (0: use htd; otherwise
                        pass to TDlib as strategy)
  -n NUMBER, --number NUMBER
                        Number of samples
  -v, --verbose         Verbose
  --model MODEL         Energy model used for sampling [bp=base pair
                        model,stack=stacking model]
  --turner              Report Turner energies of the single structures for
                        each sample
  --checkvalid          Check base pair complementarity for each structure and
                        for each sample
  --no_redundant_constraints
                        Do not add redundant constraints
  --gcweight GCWEIGHT   GC weight
  -w WEIGHT, --weight WEIGHT
                        Structure weight (def=1; in case, use last weight for
                        remaining structures)
  --plot_td             Plot tree decomposition
```

A further tool redprint_complexity.py is provided to report tree
widths for two energy models of different complexity (see our research
paper for full details) and plot the dependency graphs and tree
decompositions. This tool provides insight into the run time / space
consumption behavior of redprint on specific instances.

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

### libinfrared --- Python module exposing infrared core and extensions

Using boost::python, we expose classes of ired and ired::rnadesign to
Python, such that cluster trees can be constructed and populated with
exposed C++ constraints/functions and/or derived Python
constraints/functions.

### treedecomp.py --- Generation of tree decompositions

This python module provided by infrared supports the computation of
tree decompositions based on external libs (interfaces the Java lib
TDlib and the C++-lib libhtd; for the latter infrared implements a
wrapper (module libhtdwrap), which exposes a specific variant of
libthd-tree decomposition to Python.

### redprint.py --- Boltzmann sampling for multi-target RNA design

Redprint implements the domain-specific and higher level aspects of
Boltzmann sampling of sequences for given multi-target RNA design
instances. For either the base pair or the stacking model, it computes
the dependencies, calls the tree decomposition construction, and
constructs and populates the cluster tree for Infrared. Finally, it
calls the evaluation and sampling methods of Infrared and optionally
reports (statistics over) features of the sampled sequences.


## Known bugs / issues

* Interfacing between python and C++ can still be improved. For
  example, the boost python indexing suite is not used yet, which
  could in particular speed up Python-defined functions and
  constraints.


## State of project development --- relative stability of interfaces

The interface to the ired core is expected to remain relatively stable
by now. The domain-specific extensions in ired::rnadesign can be
expected to become more flexible in coming releases. This will cause
adaptations in the classes of redprint.py but not immediately cause
changes in their interface.
