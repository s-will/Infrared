# Example applications of Infrared

We give several examples of using Infrared to implement and solve different problems,
mostly from bioinformatics.

These examples are part of the [online documentation](https://www.lix.polytechnique.fr/~will/Software/Infrared) of Infrared and Supplemental material of the upcoming Infrared publication.

To view the notebooks, including results, visit the [application examples
in the online
documentation](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/usergroup1.html).


## How to edit and run the examples 

The examples are provided as Jupyter notebooks to support interactive
use, experimentation, and code reuse.

### Notebook file format and software requirements

The *.py files encode Jupyter notebooks as Python source in 
lightscript format. To open them as proper Jupyter notebook, install the Jupyter
and its jupytext extension. To run them, one additionally needs Infrared
and Python plotting libraries.

### Prepare software environment

We recommend the use of Conda, which allows to set up the software
environment from the command line by

```
conda create -n infrared -c conda-forge infrared jupyter jupytext matplotlib seaborn
```

Note that some notebooks, like multi-target RNA design require additional
software, which can as well be installed via Conda or using PIP.


### Start Jupyter and load notebooks

To start up the Jupyter server, activate the environment, and start jupyter

```
conda activate infrared
jupyter notebook
```

This should start the Jupyter server and open a browser window with the Jupyter file dialog; otherwise one can connect in the browser. Opening downloaded *.py files from Jupyter loads them as Jupyter notebooks.

Alternatively to using Jupyter, the notebooks can be used as Python code or from other
applications that support the Jupyter notbook format after converting them to Jupyter notebooks.
