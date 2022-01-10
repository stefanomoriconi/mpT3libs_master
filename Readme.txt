******************************************
 mpT3libs by Stefano Moriconi - Apr. 2021
******************************************

mpT3libs is a [UNIX] C-based library with a python wrapper to compute Tensors and Eigen-Decomposition in 3D for medical imaging.

The code requires OpenMP (must be installed in the computer/system) and a compiler (GCC e.g. version 6)

Compiling the shared library:
- From the folder ~/mpT3libs/, compile the shared libraries using a terminal.

In the terminal:

	cd ~/mpT3libs/
	bash compileMeFirst.sh


Once successfully compiled:

- The shared library should be in ~/mpT3libs/libmpT3libs.so
- The python wrapper is in the main directory ~/mpT3libs.py


The python wrapper consists of a set of functions useful for testing and processing:

- Tensors in 3D or to perform the Eigen-Decomposition of 3D Hessians.
- Usage and infos are in the help (comments). 

A test routine is available at ~/mpT3libs_UT/mpT3libs_test.py

Some unit tests are provided for the code in ~/mpT3libs_UT/mpT3libs_unittests.py


 -- Description --

mpT3libs package is a [UNIX] C-based library with a python wrapper to compute 3D Tensors and Eigen-Decomposition of Hessian matrices in 3D for medical imaging.

The package requires OpenMP for multi-threaded parallel computing.

The C-compiled shared library is wrapped in python for a more convenient use.

Help and documentation is given in the python wrapper (mpT3libs.py)

List of functions:

 H3_rnd: generate random Hessian
 T3_iso: generate isotropic 3D Tensor field
 T3_rnd: generate random 3D Tensor field
 msk3_to_idx: determine C-like linear indices of voxels from logical mask
 inimsk3Valid: initialise logical validity mask
 iniElvs: initialise the Eigen-Decomposition variables
 catEvs: concatenate the Eigen-Vector Components
 sepEvs: separate the Eigen-Vector Components
 rglElvs: regularise data-type for Eigen-decomposition variables
 ptrElvs: pointers to the Eigen-Decomposition variables
 ini6Cmps: initialise 6 Indep. Component variables
 ini6CmpsLE: initialise 6 Log-Euclidean Indep. Component variables
 rgl6Cmps: regularise data-type for 6 Indep. Component variables
 ptr6Cmps: pointers to the 6 Indep. Component variables
 mngIdxs: manage indexing of voxels to be processed
 mpT3_eig: (parallel) determine the Eigen-Decomposition
 mpT3_to_T3LIC: (parallel) convert 3D Tensor field to Lin. Indep. Comp.
 mpT3_to_T3LE: (parallel) convert 3D Tensor field to Log-Euclidean Comp. 
 mpT3LIC_to_T3: (parallel) convert Lin. Indep. Comp. to 3D Tensor field
 mpT3LE_to_T3: (parallel) convert Log-Euclidean Comp. to 3D Tensor field
