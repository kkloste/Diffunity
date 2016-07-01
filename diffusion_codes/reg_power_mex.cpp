/**
 * @file reg_power_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements the product of a vector with a power of a sparse matrix
 * such that the input vector can be sparse or dense, and both cases are handled
 * in a sparse manner, e.g. only nonzeros of the iterative vector are operated on
 * throughout the algorithm. The output is a dense vector.
 * The input matrix must be in matlab sparse format.
 *
 * Purdue University, June 2015.
 */



#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include <vector>
#include <assert.h>
#include <cmath>
#include <tuple>

#include "mex.h"
#include "matrix.h"

/**
 * @param G - sparse matrix values (CSC)
 * @param y - the input vector (spars,e length n)
 * @param power - number of matvecs (0 < power < inf)
 * @param npushes - the number of operations (an output)
 * @param ncols - the number of column-accesses (an output)
 * @param y - the output vector (length n)
 */ 
/* USAGE
* to compute y = A^power*v by handling iterative vectors vk in a sparse manner.
* [y npushes ncols] = reg_power_mex(P,v,power,debugflag)
*/


 
#define DEBUGPRINT(x) do { if (debugflag >= 1) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;

typedef std::pair<mwIndex,double> node_value_pair;
typedef std::vector< node_value_pair > coord_pairs;


struct sparserow {
    mwSize n, m;
    mwIndex *ai;
    mwIndex *aj;
    double *a;
    double volume;

	/**
	 * Returns the degree of node u in sparse graph s
	 */
	mwIndex sr_degree(sparserow *s, mwIndex u) {
	    return (s->ai[u+1] - s->ai[u]);
	}

};



void sparse_matvec(sparserow* G, coord_pairs& v, double* y, const mwIndex power, 
				double& totalpush, double&  totalcols)
{
    // initialize bookkeeping variables
    double total_pushes = 0.0;
    double total_columns = 0.0;
    
    // (0) Initial matvec setup 
    std::vector<mwIndex> nnz_inds(1000,0);
    mwIndex last_active = v.size();

    DEBUGPRINT(("\n             reg_power_mex: for loop starting, power = %i ", power));
    // (1) FOR LOOP for matvecs
    for (mwIndex iteration = 0; iteration < power-1; iteration++){
    DEBUGPRINT(("\n                          iteration = %i ", iteration ));
        // (2) ACTUAL MATVEC
        for( mwIndex entry = 0; entry < last_active; entry++){
            // get next nzi
            node_value_pair vpair = v[entry];
            mwIndex vrow = std::get<0>(vpair);
            double vval = std::get<1>(vpair);
        
            // perform operation
            for (mwIndex nzi=G->ai[vrow]; nzi < G->ai[vrow+1]; ++nzi) {
                mwIndex Gind = G->aj[nzi];
                double update = vval*(G->a[nzi]);
                if ( (y[Gind] == 0) && (update != 0) ){ nnz_inds.push_back(Gind); }
                y[Gind] += update;
            }
            total_columns++;
            total_pushes += (double)G->sr_degree(vrow);
        }
    DEBUGPRINT(("\n                          for loop done "));
        // (2.1) prepare next matvec
        v.clear();
        mwIndex num_nnz = nnz_inds.size();
    DEBUGPRINT(("\n                          num_nnz = %i ", num_nnz));
        for (mwIndex ind = 0; ind < num_nnz; ind++){
            mwIndex ri = nnz_inds[ind];
            double rij = y[ri];
            if ( rij != 0 ){ v.push_back( std::make_pair(ri,rij) ); }        
            y[ri] = 0;
        }
        nnz_inds.clear();
        last_active = v.size();

    DEBUGPRINT(("\n                          y scraped, v.size() = %i ", v.size() ));        
        // now y and nnz_inds are clear, ready to "receive" the next matvec.

    }// END OUTER FOR LOOP

    DEBUGPRINT(("\n             last_active = %i", last_active));
    // (3) LAST MATVEC
    for( mwIndex entry = 0; entry < last_active; entry++){
        // get next nzi
        node_value_pair vpair = v[entry];
        mwIndex vrow = std::get<0>(vpair);
        double vval = std::get<1>(vpair);
    
        // perform operation
        for (mwIndex nzi=G->ai[vrow]; nzi < G->ai[vrow+1]; ++nzi) {
            mwIndex Gind = G->aj[nzi];
            double update = vval*(G->a[nzi]);
            y[Gind] += update;
        }
        total_columns++;
        total_pushes += (double)G->sr_degree(vrow);
    }
                
    totalpush = total_pushes;
    totalcols = total_columns;
    return;
}// SPARSE_MATVEC code done


/* USAGE
* to compute y = A^power*v by handling iterative vectors vk in a sparse manner.
* [y totpush totcol] = reg_power_mex(P,v,power,debugflag)
*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 3 || nrhs > 4) {
        mexErrMsgIdAndTxt("reg_power_mex:wrongNumberArguments",
                          "reg_power_mex needs 3 to 4 arguments, not %i", nrhs);
    }
    if (nlhs > 3 || nlhs < 1) {
        mexErrMsgIdAndTxt("reg_power_mex:wrongNumberOutputs",
                          "reg_power_mex needs 1 to 3 outputs, not %i", nlhs);
    }
    if (nrhs == 4) {
        debugflag = (int)mxGetScalar(prhs[3]);
    }    
    
    DEBUGPRINT(("\n reg_power_mex: preprocessing start"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];

    // Decode sparse input matrix A
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("reg_power_mex:wrongInputMatrix",
                          "reg_power_mex needs sparse input matrix A");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("reg_power_mex:wrongInputMatrixDimensions",
                          "reg_power_mex needs square input matrix");
    }
    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);
    G.volume = (double)(G.ai[G.m]);
    
    if ( G.n != mxGetM(set) ){
        mexErrMsgIdAndTxt("reg_power_mex:wrongInputMatrixDimensions",
                          "reg_power_mex needs compatible matrix and vector dimensions");
    }    
        
    // get input vector
//    coord_pairs input_vector(1000);
    coord_pairs input_vector;
    if ( mxIsSparse(set) == true ){// Decode sparse input vector v
        sparserow v;
        v.m = mxGetM(set);
        v.n = mxGetN(set);
        v.ai = mxGetJc(set);
        v.aj = mxGetIr(set);
        v.a = mxGetPr(set);
        v.volume = (double)(v.ai[v.m]);
        for (mwIndex nzi = v.ai[0]; nzi < v.ai[1]; nzi++) {
            double rij = v.a[nzi];
            mwIndex ri = v.aj[nzi];
            if ( rij != 0 ){ input_vector.push_back( std::make_pair(ri,rij) ); }
        }
    }
    else{
        double *vi = mxGetPr(set);
        for (mwIndex ri=0; ri < G.n; ri++){
            double rij = vi[ri];
            if ( rij != 0 ){
                input_vector.push_back( std::make_pair(ri,rij) );
            }        
        }
    }    

    DEBUGPRINT(("\n            sparse data structures decoded"));
    
    double totpushes, totcols;

    mwIndex power = 1;
    if (nrhs >= 3) { power = (mwIndex)ceil(mxGetScalar(prhs[2])); }
    DEBUGPRINT(("\n            power = %i", power));    

    if ( power <= 0 ){
        mexErrMsgIdAndTxt("reg_power_mex:wrongParameterPower",
                          "reg_power_mex needs 0 < power ");
    }
        
    DEBUGPRINT(("\n            parameters obtained, everything initialized"));
    
    // dense output
    mxArray* solnvec = mxCreateDoubleMatrix(G.n,1,mxREAL);
    double *yi = mxGetPr(solnvec);
    sparse_matvec(&G, input_vector, yi, power, totpushes, totcols);
    if (nlhs > 0){ plhs[0] = solnvec; }              
 
    DEBUGPRINT(("\n             Prepare outputs.\n"));
    double* dummy = NULL;

    if (nlhs > 1){
        plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
        dummy = mxGetPr(plhs[1]);
        *dummy = totpushes;
    }
    if (nlhs > 2){
        plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
        *dummy = totcols;
        dummy = mxGetPr(plhs[2]);
    }
    DEBUGPRINT(("\n reg_power_mex:DONE   output parameters set."));
}
