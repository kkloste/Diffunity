/**
 * @file gendiff_mex.cpp
 * Implement a seeded diffusion clustering scheme.
 *
 *  Call with debugflag = 1 to display parameter values
 * before/after each call to a major function
 *
 *
 * USAGE:
 * [bestset,cond,cut,vol,y,npushes] = gendiff_mex(A,set,t,eps,debugflag)
 *
 *
 * TO COMPILE:
 *
 * if ismac
 *      mex -O -largeArrayDims gendiff_mex.cpp
 * else
 * mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims gendiff_mex.cpp
 *
 *
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif

#include "sparsehash/dense_hash_map.h"
#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include <math.h>

#include "sparsevec.hpp"
#include <unordered_set>
#include <unordered_map>
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





/*****
 *
 *          above:  DATA STRUCTURES
 *
 *          below:  CLUSTERING FUNCTIONS
 *
 ****/

/**
 *
 *  seeded_diffusion inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms; normalized so
 *                  set[i] = 1/set.size(); )
 *  output:
 *      y = exp(tP) * set
 *              with infinity-norm accuracy of eps * e^t
 *              in the degree weighted norm
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 *      max_push_count - the total number of steps to run
 *      Q - the queue data structure
 */

template <class Queue>
int seeded_diffusion(sparserow * G, sparsevec& set, sparsevec& y, const double eps,
            std::vector<double>& coeffs, const mwIndex max_push_count, Queue& Q)
{
    DEBUGPRINT(("seeded_diffusion interior: t=%f eps=%f \n", eps));
    mwIndex n = G->n;
    mwIndex N = (mwIndex)coeffs.size();
    double dummy = 1.0;
    for (mwIndex j=0; j< coeffs.size(); j++){
        dummy = dummy - coeffs[j];
        if (dummy <= eps/2.0){
            N = j-1;
            break;
        }
    }
    DEBUGPRINT(("seeded_diffusion: n=%i N=%i \n", n, N));
    
    // initialize the thresholds for pushing
    std::vector<double> psis(N,0.0);
    psis[0] = 1.0;
    for ( mwIndex j=1; j< N; j++){
        psis[j] = psis[j-1] - coeffs[j-1];
        mxAssert(psis[j] >= 0, "coefficients not properly scaled; must be nonnegative and sum to  <= 1");
    }
    
    // check eps and coeffs
    mxAssert(psis[N] <= eps/2, "coefficients input do not meet accuracy input");
    std::vector<double> pushcoeff(N,0.0);
    for (mwIndex k = 0; k < N ; k++){
        pushcoeff[k] = eps/(2*psis[k]*(double)N);
    }
    
    mwIndex ri = 0;
    mwIndex npush = 0;
    double rij = 0;
    // allocate data
    sparsevec rvec;
    
    // i is the node index, j is the "step"
    #define rentry(i,j) ((i)+(j)*n)
    
    // set the initial residual, add to the queue
    for (sparsevec::map_type::iterator it=set.map.begin(),itend=set.map.end();
         it!=itend;++it) {
        ri = it->first;
        rij = it->second;
        rvec.map[rentry(ri,0)]+=rij;
        Q.push(rentry(ri,0));
    }
    
    while (npush < max_push_count) {
        // STEP 1: pop top element off of heap
        ri = Q.front();
        Q.pop();
        // decode incides i,j
        mwIndex i = ri%n;
        mwIndex j = ri/n;
        
        double degofi = (double)G->sr_degree(i);
        rij = rvec.map[ri];
        // update yi
        y.map[i] += rij*coeffs[j];
        // update r, no need to update heap here
        rvec.map[ri] = 0;
        double update = rij/degofi;
        
        if (j == N-1) {
            // this is the terminal case, and so we add the column of A
            // directly to the solution vector y
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                //y.map[v] += update;
                y.map[v] += rij*coeffs[j];
            }
            npush += degofi;
        }
        else {
            // this is the interior case, and so we add the column of A
            // to the residual at the next time step.
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                mwIndex re = rentry(v,j+1);
                double reold = rvec.get(re);
                double renew = reold + update;
                double dv = G->sr_degree(v);
                rvec.map[re] = renew;
//                if (renew >= dv*pushcoeff[j+1] && reold < dv*pushcoeff[j+1]) {
                if (renew >= dv*pushcoeff[j] && reold < dv*pushcoeff[j]) {
                    Q.push(re);
                }
            }
            npush+=degofi;
        }
        // terminate when Q is empty, i.e. we've pushed all r(i,j) > eps*threshold[j]
        if ( Q.size() == 0) { return npush; }
    }//end 'while'
    return (npush);
}


struct greater2nd {
    template <typename P> bool operator() (const P& p1, const P& p2) {
        return p1.second > p2.second;
    }
};

void cluster_from_sweep(sparserow* G, sparsevec& p,
                        std::vector<mwIndex>& cluster, double *outcond, double* outvolume,
                        double *outcut)
{
    // now we have to do the sweep over p in sorted order by value
    typedef std::vector< std::pair<int, double> > vertex_prob_type;
    vertex_prob_type prpairs(p.map.begin(), p.map.end());
    std::sort(prpairs.begin(), prpairs.end(), greater2nd());
    
    // compute cutsize, volume, and conductance
    std::vector<double> conductance(prpairs.size());
    std::vector<mwIndex> volume(prpairs.size());
    std::vector<mwIndex> cutsize(prpairs.size());
    
    size_t i=0;
    tr1ns::unordered_map<int,size_t> rank;
    for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
         it!=itend; ++it, ++i) {
        rank[it->first] = i;
    }
    //printf("support=%i\n",prpairs.size());
    mwIndex total_degree = G->ai[G->m];
    mwIndex curcutsize = 0;
    mwIndex curvolume = 0;
    i=0;
    for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
         it!=itend; ++it, ++i) {
        mwIndex v = it->first;
        mwIndex deg = G->ai[v+1]-G->ai[v];
        mwIndex change = deg;
        for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
            mwIndex nbr = G->aj[nzi];
            if (rank.count(nbr) > 0) {
                if (rank[nbr] < rank[v]) {
                    change -= 2;
                }
            }
        }
        curcutsize += change;
        //if (curvolume + deg > target_vol) {
        //break;
        //}
        curvolume += deg;
        volume[i] = curvolume;
        cutsize[i] = curcutsize;
        if (curvolume == 0 || total_degree-curvolume==0) {
            conductance[i] = 1;
        } else {
            conductance[i] = (double)curcutsize/
            (double)std::min(curvolume,total_degree-curvolume);
        }
        //printf("%5i : cut=%6i vol=%6i prval=%8g cond=%f\n", i, curcutsize, curvolume, it->second, conductance[i]);
    }
    // we stopped the iteration when it finished, or when it hit target_vol
    size_t lastind = i;
    double mincond = std::numeric_limits<double>::max();
    size_t mincondind = 0; // set to zero so that we only add one vertex
    for (i=0; i<lastind; i++) {
        if (conductance[i] < mincond) {
            mincond = conductance[i];
            mincondind = i;
        }
    }
    //printf("mincond=%f mincondind=%i\n", mincond, mincondind);
    if (lastind == 0) {
        // add a case
        mincond = 0.0;
    }
    i = 0;
    for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
         it!=itend && i<mincondind+1; ++it, ++i) {
        cluster.push_back(it->first);
    }
    if (outcond) { *outcond = mincond; }
    if (outvolume) { *outvolume = volume[mincondind]; }
    if (outcut) { *outcut = cutsize[mincondind]; }
}

struct local_diffvec_stats {
    double conductance;
    double volume;
    double support;
    double steps;
    double eps;
    double cut;
};

/** Cluster will contain a list of all the vertices in the cluster
 * @param set - the set of starting vertices to use
 * @param coeffs - the diffusion coefficients
 * @param eps - the solution tolerance eps
 * @param p - the heatkernelpagerank vector
 * @param r - the residual vector
 * @param cluster - a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */

template <class Queue>
int hypercluster_gendiff_multiple(sparserow* G, const std::vector<mwIndex>& set,
                double eps, std::vector<double>& coeffs, sparsevec& p, sparsevec &r,
                Queue& q, std::vector<mwIndex>& cluster, local_diffvec_stats *stats)
{
    // reset data
    p.map.clear();
    r.map.clear();
    q.empty();
    DEBUGPRINT(("beginning of hypercluster \n"));
    
    size_t maxdeg = 0;
    for (size_t i=0; i<set.size(); ++i) { //populate r with indices of "set"
        assert(set[i] >= 0); assert(set[i] < G->n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = G->sr_degree(set[i]);
        r.map[set[i]] = 1./(double)(set.size()); // r is normalized to be stochastic
        //    DEBUGPRINT(("i = %i \t set[i] = %i \t setideg = %i \n", i, set[i], setideg));
        maxdeg = std::max(maxdeg, setideg);
    }
    
    DEBUGPRINT(("at last, seeded_diffusion: eps=%f \n", eps));
    
    int nsteps = seeded_diffusion(G, r, p, eps, coeffs, ceil(pow(G->n,1.5)), q);
    /**
     *      **********
     *        ***       seeded_diffusion       is called         ***
     *      **********
     */
    
    if (nsteps == 0) {
        p = r; // just copy over the residual
    }
    int support = r.map.size();
    if (stats) { stats->steps = nsteps; }
    if (stats) { stats->support = support; }
    
    // divide the probablities by their degree
    for (sparsevec::map_type::iterator it=p.map.begin(),itend=p.map.end();
         it!=itend;++it) {
        it->second *= (1.0/(double)std::max(G->sr_degree(it->first),(mwIndex)1));
    }
    
    double *outcond = NULL;
    double *outvolume = NULL;
    double *outcut = NULL;
    if (stats) { outcond = &stats->conductance; }
    if (stats) { outvolume = &stats->volume; }
    if (stats) { outcut = &stats->cut; }
    cluster_from_sweep(G, p, cluster, outcond, outvolume, outcut);
    return (0);
}

/** Grow a set of seeds via the input diffusion.
 *
 * @param G - sparserow version of input matrix A
 * @param seeds - a vector of input seeds seeds
 * @param coeffs - vector of diffusion coefficients
 * @param eps - the solution tolerance epsilon
 * @param fcond - the final conductance score of the set.
 * @param fcut - the final cut score of the set
 * @param fvol - the final volume score of the set
 */

void diffgrow(sparserow* G, std::vector<mwIndex>& seeds, double eps,
              sparsevec& p, std::vector<double>& coeffs,
            double* fcond, double* fcut, double* fvol, double* npushes)
{
    sparsevec r;
    std::queue<mwIndex> q;
    local_diffvec_stats stats;
    std::vector<mwIndex> bestclus;
    DEBUGPRINT(("gendiff_mex: call to hypercluster_heatkernel() start\n"));
    hypercluster_gendiff_multiple(G, seeds, eps, coeffs,
                                     p, r, q, bestclus, &stats);
    DEBUGPRINT(("gendiff_mex: call to hypercluster_heatkernel() DONE\n"));
    seeds = bestclus;
    *npushes = stats.steps;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
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
// [bestset,cond,cut,vol,y,npushes] = gendiff_mex(A,seed_set,coeffs,eps,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // arguments/outputs error-checking
    if (nrhs < 3 || nrhs > 5) {
        mexErrMsgIdAndTxt("gendiff_mex:wrongNumberArguments",
                          "gendiff_mex needs three to five arguments, not %i", nrhs);
    }
    if (nrhs == 5) {
        debugflag = (int)mxGetScalar(prhs[4]);
    }
    DEBUGPRINT(("gendiff_mex: preprocessing start: \n"));
    if ( nlhs > 6 ){
        mexErrMsgIdAndTxt("gendiff_mex:wrongNumberOutputs",
                          "gendiff_mex needs 0 to 6 outputs, not %i", nlhs);
    }
    
    // Retrieve Sparse Matrix input
    const mxArray* mat = prhs[0];
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("gendiff_mex:wrongInputMatrix",
                          "gendiff_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("gendiff_mex:wrongInputMatrixDimensions",
                          "gendiff_mex needs square input matrix");
    }
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    
    // Retrieve seed set input
    const mxArray* set = prhs[1];
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    
    // Prepare Coefficients
    const mxArray* coeffs_mex = prhs[2];
    mwSize ndim = mxGetNumberOfDimensions(coeffs_mex);
    const size_t *dims = mxGetDimensions(coeffs_mex);
    mxAssert(ndim <= 2, "Invalid format for input: must be k-by-1 sparse vector or k-by-2 array");
    mxAssert( std::min(dims[0], dims[1]) == 1, "Invalid format for input: coefficients must be input as a vector");
    double *vi = mxGetPr(coeffs_mex);
    mwIndex set_length = std::max(dims[0],dims[1]);
    std::vector<double> coeffs;
    for (mwIndex ind=0; ind < set_length; ind++){
        mxAssert(vi[ind] >= 0, "coefficients must be nonnegative");
        coeffs.push_back(vi[ind]);
    }

    // retrieve accuracy
    double eps = pow(10,-3);
    if (nrhs >= 4) eps = mxGetScalar(prhs[3]);

    // Prepare outputs
    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* npushes = mxCreateDoubleMatrix(1,1,mxREAL);
    
    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut; }
    if (nlhs > 3) { plhs[3] = vol; }
    if (nlhs > 5) { plhs[5] = npushes; }
    // input/outputs prepared!
    DEBUGPRINT(("gendiff_mex: preprocessing end: \n"));
    
    // execute actual code
    sparsevec diffvec;
    diffgrow(&r, seeds, eps, diffvec, coeffs, mxGetPr(cond), mxGetPr(cut), mxGetPr(vol), mxGetPr(npushes));
    
    DEBUGPRINT(("gendiff_mex: call to diffgrow() done\n"));
    
    if (nlhs > 0) { // sets output "bestset" to the set of best conductance
        mxArray* cassign = mxCreateDoubleMatrix(seeds.size(),1,mxREAL);
        plhs[0] = cassign;
        
        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<seeds.size(); ++i) {
            ci[i] = (double)(seeds[i] + 1);
        }
    }
    if (nlhs > 4) { // sets output "y" to the diffusion vector computed
        mxArray* diffvec_mex = mxCreateDoubleMatrix(r.n,1,mxREAL);
        plhs[4] = diffvec_mex;
        double *ci = mxGetPr(diffvec_mex);
        for (sparsevec::map_type::iterator it=diffvec.map.begin(),itend=diffvec.map.end();
             it!=itend;++it) {
            ci[it->first] = it->second;
        }
    }
}
