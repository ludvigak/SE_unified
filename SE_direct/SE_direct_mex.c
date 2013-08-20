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

#ifdef EWALD_REAL
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_direct_real(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "RS2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_direct_real(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "RS3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_RSRC
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_direct_real_rc(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "RC2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_direct_real_rc(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "RC3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_FD
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_direct_fd(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "FD2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_direct_fd(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG "FD3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_K0
#define EWALD_KERNEL SE2P_direct_k0(phi, idx, num_eval, x, q, N, opt)
#define EWALD_TAG " K0 "
#endif

#ifdef EWALD_SELF
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_direct_self(phi, idx, num_eval, q, N, opt)
#define EWALD_TAG "RC2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_direct_self(phi, idx, num_eval, q, N, opt)
#define EWALD_TAG "RC3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

/* common option-unpacking */
void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));
#ifdef TWO_PERIODIC
    opt->box[0] = box[0];
    opt->box[1] = box[1];
#elif THREE_PERIODIC
    opt->box[0] = box[0];
    opt->box[1] = box[1];
    opt->box[2] = box[2];
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif

    // layers: mandatory for ewald sums that are truncated 
    const mxArray* mx_layers = mxGetField(mx_opt,0,"layers");
    if(mx_layers)
	opt->layers = (int)mxGetScalar(mx_layers);
    else
	opt->layers = -1;

    // rc: mandatory for short-range real sum 
    const mxArray* mx_rc = mxGetField(mx_opt,0,"real_cutoff");
    if(mx_rc)
	opt->rc = mxGetScalar(mx_rc);
    else
	opt->rc = -1;
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
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
 
    const double* x = mxGetPr(X);
    const double* q = mxGetPr(Q);

    PHI = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
    double* restrict phi = mxGetPr(PHI);

    ewald_opts opt;
    unpack_opt(&opt, OPT);

    if(VERBOSE)
    {
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ",EWALD_TAG,N,num_eval);
	mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
		  opt.xi,opt.rc,opt.layers);
    }

    // call kernel
    EWALD_KERNEL;
    mxFree(idx);
}
