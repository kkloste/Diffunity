/**
 * @file sweepcut_mex.cpp
 *
 *
 * USAGE:
 * [bestset,conductance,cut,volume] = sweepcut_mex(A,node_array)
 *
 *
 * TO COMPILE:
 *
 * if ismac
 *      mex -O -largeArrayDims sweepcut_mex.cpp
 * else
 * mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims sweepcut_mex.cpp
 *
 *
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif

#include "sparsehash/dense_hash_map.h"
#include <vector>
#include <assert.h>
#include <math.h>

#include "sparsevec.hpp"
//#include <unordered_set>
//#include <unordered_map>
#define tr1ns std



#include <mex.h>

#define DEBUGPRINT(x) do { if (debugflag) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;

struct sparserow {
    mwSize n, m;
    mwIndex *ai;
    mwIndex *aj;
    double *a;
    
    /**
     * Returns the degree of node u in sparse graph s
     */
    mwIndex sr_degree( mwIndex u) {
        return (ai[u+1] - ai[u]);
    }
};

double get_cond( double cut, double vol, double tot_vol ){
	double temp_vol = std::min( vol, tot_vol - vol );
	if ( temp_vol == 0.0 ) return 1.0 ;
	return cut / temp_vol ;
}


/** Grow a set of seeds via the input diffusion.
 *
 * @param G - sparserow version of input matrix A
 */

void sweepcut(sparserow* G, std::vector<mwIndex>& noderank, mwIndex& bindex,
	double& bcond, double& bvol, double& bcut)
{
	mwIndex length_array = noderank.size();
	double total_volume = (double)G->ai[G->m] ;
    sparsevec local_cut_G;
	double  curcond = 0.0, curvol = 0.0, curcut = 0.0;
	bcond = 2.0;

	for ( mwIndex counter = 0, curindex = 0; counter < length_array; counter++, curindex++ ){
		mwIndex curinterior = 0;
		curnode = noderank[curindex];
		// move next neighbor into the local cut graph
		// and update current volume and current cut-edges
		for (mwIndex nzi = G->ai[curnode]; nzi < G->ai[curnode+1]; nzi++){
			mwIndex curneighb = G->aj[nzi];
			if (local_cut_G.map[curneigb] == 1){
				curinterior += 1;
			}
			local_cut_G.map[curneighb] = 1;
		}
		curvol += (double) G->sr_degree(curnode);
		curcut += (double) ( G->sr_degree(curnode) - 2*curinterior );
		curcond = get_cond( curcut, curvol, total_volume);
		if (curcond < bcond) {
			bcond = curcond;
			bcut = curcut;
			bvol = curvol;
			bindex = curindex;
		}
		if ( curvol >= total_volume/2.0 ){ break; }
	}
}

void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    
    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex)elem - 1;
    }
}


// USAGE
// [bestset,cond,cut,vol,y,npushes] = sweepcut_mex(A,seed_set,coeffs,eps,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // arguments/outputs error-checking
    if ( nrhs > 3 || nrhs < 2 ) {
        mexErrMsgIdAndTxt("sweepcut_mex:wrongNumberArguments",
                          "sweepcut_mex needs two to three arguments, not %i", nrhs);
    }
    if (nrhs == 3) {
        debugflag = (int)mxGetScalar(prhs[2]);
    }
    DEBUGPRINT(("sweepcut_mex: preprocessing start: \n"));
    if ( nlhs > 4 ){
        mexErrMsgIdAndTxt("sweepcut_mex:wrongNumberOutputs",
                          "sweepcut_mex needs 0 to 4 outputs, not %i", nlhs);
    }
    
    // Retrieve Sparse Matrix input
    const mxArray* mat = prhs[0];
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("sweepcut_mex:wrongInputMatrix",
                          "sweepcut_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("sweepcut_mex:wrongInputMatrixDimensions",
                          "sweepcut_mex needs square input matrix");
    }
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    
    // Retrieve ordered node array
    const mxArray* set = prhs[1];
    std::vector< mwIndex > node_array;
    copy_array_to_index_vector( set, node_array );
    
// INPUTS DECODED
// now perform sweep procedure

    // Prepare outputs
    mxArray* conductanceMX = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cutMX = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* volumeMX = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* break_indexMX = mxCreateDoubleMatrix(1,1,mxREAL);
	double conductance, cut, volume;
	mwIndex break_index = 0;
    if (nlhs > 0) { plhs[0] = break_indexMX; }
    if (nlhs > 1) { plhs[1] = conductanceMX; }
    if (nlhs > 2) { plhs[2] = cutMX; }
    if (nlhs > 3) { plhs[3] = volumeMX; }

    // input/outputs prepared!
    DEBUGPRINT(("sweepcut_mex: preprocessing end: \n"));
    
    // execute actual code
    sweepcut(&r, node_array, break_index, conductance, volume, cut);
    
    DEBUGPRINT(("sweepcut_mex: call to sweepcut() done\n"));
    
	*break_indexMX = (double)break_index + 1.0; // shifting to match Matlab indexing
	*conductanceMX = conductance;
	*cutMX = cut;
	*volumeMX = volume;
}
