Pyomo package for probabilistic design space characterisation. 

Installation instructions:
1. Install prerequisistes
    conda install numpy pyomo openmpi mpi4py
    pip install mpi-sppy
2. Build the package source:
    python setup.py sdist
3. Install the package
    pip install ../dist/pyds-.tar.gz
4. (Optional) Clean the installation directory
    python setup.py clean


