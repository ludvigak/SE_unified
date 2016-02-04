#include "mex.h"
#include "SE_Stokes_direct.h"

#define IDX prhs[0] // Indices in x to evaluate at
#define X   prhs[1] // Source locations
#define F   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define U plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

// Select what to compute from command-line defined variables.
// Typical compiler defines: 
// (...) -DEWALD_REAL -DTHREE_PERIODIC -DHASIMOTO
#ifdef EWALD_REAL
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_Stokes_direct_real(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "RS2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_Stokes_direct_real(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "RS3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_RSRC
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_Stokes_direct_real_rc(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "RC2P"
#elif THREE_PERIODIC
#ifdef EXTERNAL
#define EWALD_KERNEL SE3P_Stokes_direct_real_ext_rc(u, xt, num_eval, x, f, N, opt)
#define EWALD_TAG "RC3P_EXT"
#define XT prhs[0] // Target locations
#else
#define EWALD_KERNEL SE3P_Stokes_direct_real_rc(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "RC3P"
#endif
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_FD
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_Stokes_direct_fd(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "FD2P"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_Stokes_direct_fd(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG "FD3P"
#else
#error "Must give -D<TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_K0
#define EWALD_KERNEL SE2P_Stokes_direct_k0(u, idx, num_eval, x, f, N, opt)
#define EWALD_TAG " K0 "
#endif

#ifdef EWALD_SELF
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_Stokes_direct_self(u, idx, num_eval, f, N, opt)
#define EWALD_TAG "SELF"
#elif THREE_PERIODIC
#define EWALD_KERNEL SE3P_Stokes_direct_self(u, idx, num_eval, f, N, opt)
#define EWALD_TAG "SELF"
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

    // rc: mandatory for RSRC ewald sum 
    const mxArray* mx_rc = mxGetField(mx_opt,0,"rc");
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
    // input dims, assumes X and F is N-by-3
    const int N = mxGetM(X);
    
#ifdef EXTERNAL
    const int num_eval = mxGetM(XT);
    const double* xt = mxGetPr(XT);
#else
    const int num_eval = mxGetN(IDX); // FIXME: indices assumed to be row vec 
    const double* idx_d = mxGetPr(IDX);
    int* idx = mxMalloc(num_eval*sizeof(int)); // temporary index vector
    index_translation(idx, idx_d, num_eval);
#endif


#ifndef EWALD_SELF
    const double* x = mxGetPr(X);
#endif
    const double* f = mxGetPr(F);

    U = mxCreateDoubleMatrix(num_eval, 3, mxREAL); // the result
    double* restrict u = mxGetPr(U);

    ewald_opts opt; // parameters
    unpack_opt(&opt, OPT);

    if(VERBOSE)
    {
	mexPrintf("[STOKESLET EWALD (%s)] MEX N=(%d,%d) ",EWALD_TAG,N,num_eval);
	mexPrintf("xi = %.2f layers=%d\n",opt.xi,opt.layers);
    }

    // call kernel
    EWALD_KERNEL;
#ifndef EXTERNAL
    mxFree(idx);
#endif 
}
