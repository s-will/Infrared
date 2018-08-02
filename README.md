# Infrared

Infrared is a generic C++/Python hybrid library for efficient (fixed-parameter tractable) Boltzmann sampling. 

## Disclaimer and license
Infrared is free software. Note that the system is still in an early stage of development and is likely to undergo significant changes. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details [<http://www.gnu.org/licenses/>].

## Main features of Infrared
Infrared provides a fast and flexible C++ engine that evaluates a constraint network consisting of variables, multi-ary functions, and multi-ary constraints. Functions and constraints are C++ or Python objects, where
new functions and contraints are easily added in C++ or in Python! The evaluation is performed efficiently using cluster tree elimination following a (hyper-)tree decomposition of the dependencies (due to functions and constraints).
Interpreting the evaluations as partition functions, the system supports sampling of variable assignments from
the corresponding Boltzmann distribution.

While the core library is agnostic to the origin of the tree decomposition, it comes with interfaces to TDlib and libhtd.

## RNARedPrint 2
On top of the Infrared library, we re-implemented RNARedPrint [<https://github.com/yannponty/RNARedPrint>], which is shipped together with the library and serves as first non-trivial example for the use of the Infrared engine. RNARedprint samples RNA sequences from a Boltzmann distribution based on the energies of multiple target structures and GC-content. Our research paper on RNARedPrint "Stefan Hammer, Yann Ponty, Wei Wang, Sebastian Will.  Fixed-Parameter Tractable Sampling for RNA Design with Multiple Target Structures. Proc. of RECOMB, 2018." decribes many of the ideas that lead to the development of Infrared and points to the origins in cluster tree decomposition and constraint networks. If you find this software useful for your own work, please do not forget to cite us.

## Installation
The software can be compiled and installed (after cloning the repository) on GNU/Linux or other Unix-like systems as follows.
In preparation, install libhtd <https://github.com/mabseher/htd/releases>. Moreover, compilation requires boost.graph and boost.python. Compile and install Infrared itself by
```
cd Infrared
./autoreconf -i
./configure --with-htd=$path_to_htd_installation --prefix=$path_to_infrared_installation
make && make install
```
To use tree decomposition based on the Java tool TDlib, please obtain a copy and set the Java CLASSPATH accordingly.

## Running Redprint
For running RNARedprint, one furthermore needs an installation of the ViennaRNA package; make sure that the Python module is found in the path described by the environment variable PYTHONPATH. PYTHONPATH must as well point to the libinfrared module from the Infrared installation, i.e. $path_to_infrared_installation/lib

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

Generate training data

positional arguments:
  infile                Input file

optional arguments:
  -h, --help            show this help message and exit
  --method METHOD       Method for tree decomposition (0: use htd; othwerwise
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
A further tool redprint_complexity.py is provided to report tree widths for two energy models of different complexity (see our research paper for full details) and plot the dependency graphs and tree decompositions. This tool provides insight into the run time / space consumption behavior of redprint on specific instances.

## Known bugs / issues
* By default redprint.py performs tree decomposition using libhtd, however we do not make optimal use of the library. In
  result, the efficiency can strongly vary with the non-deterministically generated and non-optimized tree decomposition.
* Currently, the engine still requires redundant constraint insertions for improving its performance
  (as demonstrated by the redprint implementation). This should be avoidable by suitable optimizations.
* Interfacing between python and C++ can still be improved. For example, the boost python indexing suite is not used yet,
  which could in particular speed up Python-defined functions and constraints.
