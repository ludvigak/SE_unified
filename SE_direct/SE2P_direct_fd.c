#include <math.h>
#include "mex.h"

#define IDX prhs[0]
#define X   prhs[1] // Source locations
#define Q   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define PHI plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define __EXP_ARG_MAX 600
#define PI 3.141592653589793

typedef struct 
{
    double box[3];
    double xi;
    int layers; 
    double rc; 
} ewald_opts;

void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));

    opt->box[0] = box[0];
    opt->box[1] = box[1];

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

static inline double theta_plus(double z, double k, double xi)
{
    /* idea for a more stable form [LK] */
    /* exp( k*z + log( erfc(k/(2.0*xi) + xi*z) ) ); */

    if(k*z <  __EXP_ARG_MAX)
	return exp( k*z)*erfc(k/(2.0*xi) + xi*z);
    else 
	return 0.0;
}

static inline double theta_minus(double z, double k, double xi)
{
    /* exp(-k*z + log( erfc(k/(2.0*xi) - xi*z) ) ); */

    if(-k*z <  __EXP_ARG_MAX)
	return exp(-k*z)*erfc(k/(2.0*xi) - xi*z);
    else 
	return 0.0;
}

void SE2P_direct_fd(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
    double k[2], xm[3]; 
    double kn, k_dot_r, z, phi_m;
    double cm, cp;
    const double xi = opt.xi;

    for(int m=0; m<nidx; m++)
    {
	xm[0] = x[idx[m]    ];
	xm[1] = x[idx[m]+N  ];
	xm[2] = x[idx[m]+2*N];
	phi_m = 0;
	for(int n = 0; n<N; n++){
	    for(int j0 = -opt.layers; j0<=opt.layers; j0++)
		for(int j1 = -opt.layers; j1<=opt.layers; j1++)
		{
		    if(j0 == 0 && j1 == 0)
			continue;

		    k[0] = 2*PI*j0/opt.box[0];
		    k[1] = 2*PI*j1/opt.box[1];
		    kn = sqrt(k[0]*k[0] + k[1]*k[1]);
		    k_dot_r = k[0]*(xm[0]-x[n]) + k[1]*(xm[1]-x[n+N]);
		    z = xm[2]-x[n+2*N];
		    cp = theta_plus(z,kn,xi);
		    cm = theta_minus(z,kn,xi);

		    phi_m += q[n]*cos(k_dot_r)*(cm+cp)/kn;
		}
	}
	phi[m] = PI*phi_m/(opt.box[0]*opt.box[1]);
    }
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
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","FD2P",N,num_eval);
	mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
		  opt.xi,opt.rc,opt.layers);
    }

    // call kernel
    SE2P_direct_fd(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}
