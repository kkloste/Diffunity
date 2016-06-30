# Generalized Diffusion method


Description of the files in the `/code/` directory.


## `/diffusion_codes`

Contains C, MEX, and matlab codes for computing diffusions, including our general diffusion push algorithm, the ACL ppr push algorithm, and a fast sparse method for computing P^k * v. To use any of these codes, you must first call:

* `compile.m` -- call this function from `[project]/code/diffusion_codes/` to compile all MEX codes.

Codes for algorithms:

* `pprgrow_mex.cc` -- Fast implementation of ACL push for ppr, from [Gleich & Seshadrhi 12].
* `pprgrow.m` -- matlab wrapper to call `pprgrow_mex.cc`
* `gendiff_mex.cpp` -- MEX code for our novel push algorithm for generalized diffusions.
* `gendiff_grow.m` -- NOT DONE: will be a matlab wrapper for gendiff

## `/test`

Contains Matlab test scripts.

* `test_gendiff_mex.m` -- after compiling our diffusion codes, can check that `gendiff_mex.cpp` is working by calling this script. This computes PageRank explicitly, and also using `pprgrow_mex` and `gendiff_mex`, and outputs the error.


## `/ls`

### Experiment codes
1. In the directory `/experiments/comm_analysis`, files `community_analysis`, `comm_tprs`, and `comm_analysis2` compute a number of structural statistics for all the communities of all the datasets.


### Code for diffusions:
1. `gendiff_mex.cpp` implements our generalized push method.
1. `hpolymv.m` computes a vector times a polynomial of a matrix. This is intended specifically to serve as true solution vectors to compare against our new method.
2. `pprgrow_mex.cc` and `pprgrow.m` are from [ Gleich & Seshadhri KDD 2012 ]
	* adapted to output both the personalized pagerank diffusion and its set of best conductance.
	
### Other code:
* `test_gendiff_mex.m` confirms that gendiff_mex outputs a solution vector satisfying the desired accuracy requirement. It also compares the accuracy attained by gendiff with the accuracy attained by running pprgrow with the same parameter settings.
