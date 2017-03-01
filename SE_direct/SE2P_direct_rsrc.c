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

void SE2P_direct_real_rc(double* restrict phi, 
			 const int* restrict idx, int nidx,
			 const double* restrict x, 
			 const double* restrict q, int N,
			 const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
    for(int m=0; m<nidx; m++)
    {
	p = 0;
	for(int n=0; n<N; n++)
	{
	    rvec[0] = x[idx[m]    ]-x[n    ];
	    rvec[1] = x[idx[m]+N  ]-x[n+N  ];
	    rvec[2] = x[idx[m]+2*N]-x[n+2*N];
	    qn = q[n];

	    for(int p0 = -opt.layers; p0<=opt.layers; p0++)
		for(int p1 = -opt.layers; p1<=opt.layers; p1++)
		{
		    if(idx[m] == n && p1 == 0 && p0 == 0)
			continue;
			
		    r = sqrt((rvec[0]+p0*opt.box[0])*
			     (rvec[0]+p0*opt.box[0])+
			     (rvec[1]+p1*opt.box[1])*
			     (rvec[1]+p1*opt.box[1])+
			     rvec[2]*rvec[2]);
			
		    if(r > opt.rc) continue;

		    p += qn*erfc(opt.xi*r)/r;
		}
	}
	phi[m] = p;
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
      mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","RSRC2P",N,num_eval);
      mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
		opt.xi,opt.rc,opt.layers);
    }

    // call kernel
    SE2P_direct_real_rc(phi, idx, num_eval, x, q, N, opt);

    mxFree(idx);
}
