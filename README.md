# WMM
Wavefront Marching Methods: a unified algorithm to solve Eikonal and static Hamilton-Jacobi equations

## Instructions for compiling:
Run the following commands:

    cd python/cython_files/two_dim
    python setup.py build_ext --inplace 
    cd ../three_dim
    python setup.py build_ext --inplace 
    cd ../distance
    python setup.py build_ext --inplace
