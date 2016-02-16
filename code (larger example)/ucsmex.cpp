#include "mex.h"
#include "sparsifier.h"
#include <cassert>

// compile with
// mex ucsmex.cpp sparsifier.cpp solver.cpp -lmwblas -lmwlapack
// from matlab terminal

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[] )
{
// prhs - array of input arguments
// plhs - array of output arguments
// nrhs - number of input arguments
// nlhs - number of output arguments
//
// expected inputs:
//      U - data matrix
//      k - problem rank
//      l - oversampling size
//
    
    if (nrhs < 3) { mexPrintf("Insufficient data inputs\n"); return; }
    
    // input data
    double* U;
    U = mxGetPr(prhs[0]);
    int n = mxGetM(prhs[0]); // gets number of rows of a matrix
    int k = mxGetScalar(prhs[1]);
    int l = mxGetScalar(prhs[2]);
    assert (k<l && l<=n);

    // run column selection
    sparsifier S(n,k,l);
    double* pi;
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    pi = mxGetPr(plhs[0]);
    for (int i=0; i<n; ++i) { pi[i] = 0; }
    S.columnSelect(pi, U);
    return;
}


