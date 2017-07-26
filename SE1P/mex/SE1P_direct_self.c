#include "mex.h"
#include "SE_direct.h"

#define IDX prhs[0]
#define X   prhs[1] // Source locations
#define Q   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define PHI plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

/* common option-unpacking */
void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
}

void SE1P_direct_self(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict q, int N, 
		      const ewald_opts opt)
{
    double c = 2*opt.xi/sqrt(PI);
    for(int m=0; m<nidx; m++)
	phi[m] = -c*q[idx[m]];
}

/* no input checking is done */
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    // input dims
    const int N = mxGetM(X);
    
    const int num_eval = mxGetN(IDX); // FIXME: indices assumed to be row vec
    const double* idx_d = mxGetPr(IDX);
    int* idx = mxMalloc(num_eval*sizeof(int));
    index_translation(idx, idx_d, num_eval);
 
    const double* q = mxGetPr(Q);

#ifndef FORCE
    PHI = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
    double* restrict phi = mxGetPr(PHI);
#else 
    /* This is to allocate 3 vectors for the force. 
     * (FIXME) Note that the variable is still called PHI.*/
    PHI = mxCreateDoubleMatrix(num_eval, 3, mxREAL);
    double* restrict phi = mxGetPr(PHI);
#endif

    ewald_opts opt;
    unpack_opt(&opt, OPT);

    if(VERBOSE)
    {
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","SELF1P",N,num_eval);
	mexPrintf("xi = %.2f\n",opt.xi);
    }

    // call kernel
    SE1P_direct_self(phi, idx, num_eval, q, N, opt);
    mxFree(idx);
}

