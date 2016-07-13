# Diffunity -- diffusion codes and related utilities with applications to community detection


## TO DO

* Implement the Andersen Lang truncated lazy random walk algorithm.
* write special version of `gendiff_mex` designed specifically to compute P^k*s .



## `/diffusion_codes`

Contains C, MEX, and matlab codes for computing diffusions, including our general diffusion push algorithm, the ACL ppr push algorithm, a heat kernel push algorithm, and a fast sparse method for computing P^k * v. To use any of these codes, you must first call:

* `compile.m` -- call this function from `[project]/diffusion_codes/` to compile all MEX codes.

Codes for algorithms:

1. `gendiff_mex.cpp` implements our generalized push method.
* `gendiff_grow.m` -- NOT DONE: will be a matlab wrapper for gendiff
2. `pprgrow_mex.cc` and `pprgrow.m` are from [ Gleich & Seshadhri KDD 2012 ]
	* adapted to output both the personalized pagerank diffusion and its set of best conductance.
* `reg_power_mex.cpp` -- given sparse input matrix A and vector v, computes A^k*v exactly (i.e. no rounding) but in a sparse manner for an input power k.
* `hkgrow_mex.cpp` -- up-to-date implementation of heat kernel diffusion from [Gleich & Kloster 2014]


* `sweepcut_mex.cpp` -- fast algorithm for sweeping over an input vector.

## `/test`

Contains Matlab test scripts.

* `test_gendiff_mex.m` -- after compiling our diffusion codes, can check that `gendiff_mex.cpp` is working by calling this script. This computes PageRank explicitly, and also using `pprgrow_mex` and `gendiff_mex`, and outputs the error.


### Experiment codes





### Other code:
* `test_gendiff_mex.m` confirms that gendiff_mex outputs a solution vector satisfying the desired accuracy requirement. It also compares the accuracy attained by gendiff with the accuracy attained by running pprgrow with the same parameter settings. confirms that gendiff_mex outputs a solution vector satisfying the desired accuracy requirement. It also compares the accuracy attained by gendiff with the accuracy attained by running pprgrow with the same parameter settings.ow with the same parameter settings.

* `test_sweepcut.m` -- to check that `sweepcut_mex` is operating correctly.
