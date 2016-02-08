#include "mex.h"
#include "rotlet_direct.h"

#define XT  prhs[0] // Target locations
#define X   prhs[1] // Source locations
#define F   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define U plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

/* common option-unpacking */
void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));
    opt->box[0] = box[0];
    opt->box[1] = box[1];
    opt->box[2] = box[2];
    const mxArray* mx_rc = mxGetField(mx_opt,0,"rc");
    mxAssert(mx_rc > 0, "rc must be positice");
    opt->rc = mxGetScalar(mx_rc);
}


void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int Nt = mxGetM(XT);
    const double* xt = mxGetPr(XT);
    const int N = mxGetM(X);
    const double* x = mxGetPr(X);
    const double* f = mxGetPr(F);

    U = mxCreateDoubleMatrix(Nt, 3, mxREAL); // the result
    double* restrict u = mxGetPr(U);

    ewald_opts opt; // parameters
    unpack_opt(&opt, OPT);

    // call kernel
    rotlet_direct_rsrc(u, xt, Nt, x, f, N, opt);
}
