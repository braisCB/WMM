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
   
## Citing:
If you use this code, please cite:

    @article{cancela2020wavefront,
        title={Wavefront Marching Methods: a unified algorithm to solve Eikonal and static Hamilton-Jacobi equations},
        author={Cancela, Brais and Alonso-Betanzos, Amparo},
        journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
        year={2020},
        publisher={IEEE}
    }
