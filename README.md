[TOC]

# Infrared
![](Doc/Infrared-Logos/infrared-small.gif)

Infrared is a framework for efficient, tree decomposition based solving of
declaratively modeled problems. Models can be solved by optimization or
Boltzmann sampling. The latter allows targeting of features by
multi-dimensional Boltzmann sampling and supports further stochastic
optimization.

## Resources

**Online documentation**: <https://www.lix.polytechnique.fr/~will/Software/Infrared>

**Usage examples and tutorials**: Various examples of solving bioinformatics problems in
Infrared provide an important entry point for practical use of the system
(see online Documentation, Examples). Tutorials are provided for RNA design
in Infrared. One is accompanying a book chapter and design in Infrared; the
other, directly showcasing elementary usage.

**Preprint**: [HT Yao, B Marchand, SJ Berkemer, Y Ponty, S Will. 
Infrared: a declarative tree decomposition-powered framework for bioinformatics. 2023](https://inria.hal.science/hal-04211173)

This manuscript comprehensively describes the system and also discusses many of the given example applications.

**Installation**: [```conda install -c conda-forge infrared```](https://anaconda.org/conda-forge/infrared)

<table>
<thead>
<tr>
<th>Name</th>
<th>Downloads</th>
<th>Version</th>
<th>Platforms</th>
</tr>
</thead>
<tbody>
<tr>
<td><a href="https://anaconda.org/conda-forge/infrared" rel="nofollow"><img src="https://camo.githubusercontent.com/c930002e4cb49ed5db671880df00436fa9938b06652a510f3e9c39eb33f5461e/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f7265636970652d696e6672617265642d677265656e2e737667" alt="Conda Recipe" data-canonical-src="https://img.shields.io/badge/recipe-infrared-green.svg" style="max-width: 100%;"></a></td>
<td><a href="https://anaconda.org/conda-forge/infrared" rel="nofollow"><img src="https://camo.githubusercontent.com/f3a395ec56108843b3283a11221fe7d4c21967acf623cc31af195f3d50625b68/68747470733a2f2f696d672e736869656c64732e696f2f636f6e64612f646e2f636f6e64612d666f7267652f696e6672617265642e737667" alt="Conda Downloads" data-canonical-src="https://img.shields.io/conda/dn/conda-forge/infrared.svg" style="max-width: 100%;"></a></td>
<td><a href="https://anaconda.org/conda-forge/infrared" rel="nofollow"><img src="https://camo.githubusercontent.com/66b47863879e476698c41b7ab784b5f018d5db57a05a89201b3b359e1bc7a295/68747470733a2f2f696d672e736869656c64732e696f2f636f6e64612f766e2f636f6e64612d666f7267652f696e6672617265642e737667" alt="Conda Version" data-canonical-src="https://img.shields.io/conda/vn/conda-forge/infrared.svg" style="max-width: 100%;"></a></td>
<td><a href="https://anaconda.org/conda-forge/infrared" rel="nofollow"><img src="https://camo.githubusercontent.com/4fc990acda82ad8788d0b8ae02468d231b7d4fc741947314a25a8a232189266d/68747470733a2f2f696d672e736869656c64732e696f2f636f6e64612f706e2f636f6e64612d666f7267652f696e6672617265642e737667" alt="Conda Platforms" data-canonical-src="https://img.shields.io/conda/pn/conda-forge/infrared.svg" style="max-width: 100%;"></a></td>
</tr>
</tbody>
</table>

**Gitlab repository**: <https://gitlab.inria.fr/amibio/Infrared>


## Main features

Infrared solves a broad class of sampling and optimization problems that
can be expressed as feature networks by efficient tree decomposition based
algorithms. Solving is thus performed with parameterized complexity
depending on the treewidth of the network.  The system was developed to
target bioinformatics problems with potentially complex dependencies,
originally RNA sequence design with multiple target structures.  Such
problems can be specified in a declarative, compositional style as
(evaluated) constraint models using the Python high-level interface of
Infrared.

Accessible through the Python interface, the framework provides a fast and
flexible C++ engine that evaluates constraint networks consisting of
variables, multi-ary functions, and multi-ary constraints.  This evaluation
by generic algorithms is efficient depending on the complexity, measured as
tree-width, of the network. 

Application specific functions and constraints can be directly defined in
Python.  The evaluation is performed efficiently using cluster tree
elimination following a (hyper-)tree decomposition of the dependencies (due
to functions and constraints). Interpreting the evaluations as partition
functions, the system supports sampling of variable assignments from the
corresponding Boltzmann distribution.  Finally, Infrared implements a
generic multi-dimensional Boltzmann sampling strategy to target specific
feature values. Such functionality is made conveniently available via
general Python classes. In particular, the interface was defined to allow
the straightforward and declarative specification of the feature network
model.

## Installation

### Conda installation

Infrared is installed most easily using conda. Infrared is depolyed on conda-forge channel.
Users can skip the first line command if it's already done.

```
conda config --add channels conda-forge 
conda install infrared
```

### Pip installation from source 

For users who don't want to use conda, Infrared can also be installed with standard pip install from
it's source, which we make
freely available in [Infrared's Gitlab repository](https://gitlab.inria.fr/amibio/Infrared).
Compiling and installing requires a C++ / Python
build environment including cmake, and installation of further dependencies, e.g. [pybind11](https://github.com/pybind/pybind11)
and [Treedecomp](https://gitlab.inria.fr/amibio/treedecomp).

After installing dependencies, one compiles and installs Infrared from its base directory by
```
python3 -m pip install .
```

Treedecomp can be installed analogously; Infrared requires at least version 1.1.0.

Note that older Linux distributions, e.g. Ubuntu 18.04, install only outdated versions of pybind11 via their package managers; we require at least version 2.4. One can install pybind11 as well via pip by 
```
PYBIND11_GLOBAL_SDIST=1 python3 -m pip install https://github.com/pybind/pybind11/archive/master.zip
```

## Documentation

We provide [Infrared's documentation](https://www.lix.polytechnique.fr/~will/Software/Infrared) online. The documentation comprises general information, API reference and examples for the use of Infrared's high-level Python interface.

Jupyter notebooks with code examples are part of the online documentation and can be as well downloaded from subdirectory [Doc](https://gitlab.inria.fr/amibio/Infrared/-/tree/master/Doc) on Inrared's Gitlab repository.

A further entry-point to using the library in novel sampling applications
is provided by the code of [RNARedPrint 2](https://gitlab.inria.fr/amibio/RNARedPrint)
and [RNAPOND](https://gitlab.inria.fr/amibio/RNAPOND).

## Infrared architecture and background

The system was originally build to separate the core "Infrared" from applications
like [RNARedprint 2](https://gitlab.inria.fr/amibio/RNARedPrint) or
[RNAPOND](https://gitlab.inria.fr/amibio/RNAPOND).


### Python high level interface

We expose classes of the C++ library to Python, such that problems can be
implemented (i.e. modeled) and solved by writing Python code.  Python
programs  using the infrared library, typically import module infrared,
create an instance of Model and populate it with variables (specifying
their finite domains), constraints and functions. The model automatically
generates features from the functions; in order to generate and control
several features, the functions can be assigned to function groups.
Moreover, the user can define additional features.

The model is passed to Solvers, i.e. Optimizers or Samplers that generate samples
from a specific Boltzmann distrubition as well as samples targeted at
specific feature values. The latter performed by multi-dimensional
Boltzmann sampling. The module provides access to the Infrared core and
exports the additional base classes 

Internally, as well on the Python side, cluster trees are constructed and
populated with constraints and functions. 

### Infrared C++ core

The C++ component solves feature networks, i.e. optimizes or samples, based
on cluster trees.

A /cluster tree/ of a feature network (of variables,
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

The core preforms precomputation, i.e. calculation of partition
functions, and Boltzmann sampling based on a given, populated cluster
tree.

By design, the Infrared core is strictly low-level and
domain-agnostic, it only knows finite domain variables, constraints,
and functions; as well as to evaluate partition functions and sample
based on cluster trees for constraint networks. This will allow using
the same core for evaluation (even not necessarily of partition
functions, e.g. instead optimizing energy/fitness) and sampling in
very different domains just by writing domain-specific Python code.


## Disclaimer and license

Infrared is free software. It was part of project [RNARedPrint](https://github.com/yannponty/RNARedPrint), then separated as a stand alone project.
Note that the system is in active development and is likely to still undergo 
changes. It is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details
[<http://www.gnu.org/licenses/>].

