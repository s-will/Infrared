ChangeLog
=========

Release 1.0b - 2021-11-16
--------------------------

This release facilitates the modeling of sampling and optimization problems 
from a much broader range of applications by

* introducing a new modeling syntax which facilitates the definition
  of constraint models and supporting the definition of 
  new constraints and functions with simple Syntax in Python

* enabling modeling in a truly declarative and compositional style

* supporting domain bounds, which can
  be performance critical for many applications

* optimizing space requirements

* enabling optimization in addition to sampling, where
  both can be based on the same constraint model

* improving the API documentation

* adding tutorials and examples in Jupyter notebooks that
  demonstrate the modeling of design problems in the framework

### Python interface

* Introduce new modeling syntax to define constraint network models based on a
  new class ```Model```

* Old class ConstraintNetwork is merged with new class Model

* Support targeting of only selected features

* Revise syntax for the definition of function and constraints

* Add support for optimization (in addition to sampling)
  - allow to select evaluation algebra
  - introduce weighted functions

* Support IUPAC sequence constraints with a specific form of
  constraint simplification

* Rename functions and constraints in module rna

* Support named variables in constraint models

### C++ library

* Support upper and lower bound of variable domains; add class
  ```FiniteDomain```

* Remove dependency on ```boost::graph```

* Implement traceback for constraint optimization

* Optimize space consumption of materialized functions and
constraints by using a specialized replacement
```SimpleMap``` of standard sparse containers like
```std::hash_map```

* Add evaluation in arctic semiring for constraint optimization

### Documentation

* Revisit API documentation

* Adapt tutorial notebook to new Python interface

* Add example for multi target design from Python code (Jupyter
  notebook)

* Add Jupyter notebook as supplement of a book chapter on Infrared
  with new example code

### Bug fixes

* Fix bug of repeated dependencies in the specification of models

* Fix undefined behavior in the materialization of
functions/constraints


Release 0.5.1 - 2021-07-10 
--------------------------

### Python Interface

* Export algebra for constraint optimization

* Redefine constraints and functions in Python

### C++ library

Materialization of functions /and/ constraints

* supports the definition of functions and constraints in
  Python without performance penalties

* the existing mechanism for functions was extended to constraints

### Bug Fix

* Update behavior of ```BoltzmannSampler.plot_td()```

### Documentation

* Re-enable API documentation via doxygen

* Update tutorial jupyter notebook

Release 0.5 - 2021-05-04
------------------------

### C++ library
Use ```pybind11```

### Installation and dependencies

* Change build mechanism to ```setup.py``` and ```cmake```

* Enable installation via pip and conda

Release 0.4 - 2020-11-08
------------------------

### Python Interface

* Refactor and streamline the Python interface

* Support tree decomposition by networkx, improve the tree decomposition interface

* Add access to partition function and provide consistency check in

### Documentation

* Add tutorial of the Python interface (Jupyter notebook)

### C++ library

* Add tests, fix handling of inconsistent networks

* Revert to using ```boost::python``` (to resolve a bug in the use of ```pybind11```)

### Installation and dependencies

* Improve installation via autotools

* Replace dependency on ```libhtd``` by ```networkx```

Release 0.3 - 2018-08-17
------------------------

* Generalize Python classes for Boltzmann sampler construction and
  provide them in Python module ```infrared``` for simpler use of the lib

* Implement multi-dimensional Boltzmann sampling and provide this
  functionality via module ```infrared```

Release 0.2.1 - 2018-08-08
--------------------------

* Add/cleanup Doxygen documentation of Python modules

Release 0.2 - 2018-08-07
------------------------

* Improve program comments, support doxygen (```make doxygen-doc```)

* Revise code, fix make check

Release 0.1 - 2018-08-07
------------------------

* Optimize infrared evaluation/sampling engine

* Optimize tree decomposition using ```libhtd```

* Cleanup treedecomp and redprint Python code and add features

* Bugfixes

* Implement some tests

Release 0.0.1 - 2018-08-02
--------------------------

Initial release
