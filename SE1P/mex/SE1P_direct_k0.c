#include "mex.h"
#include "SE_direct.h"
#include "mathint.h"

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
    if(opt->xi==0)
      mexErrMsgTxt("xi cannot be zero");
    
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));

    opt->box[2] = box[2];
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
}

#ifdef FORCE
void SE1P_direct_k0(double* restrict force,
		    const int* restrict idx, int nidx,
		    const double* restrict x,
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
  double rho2,xm[2];
  const double xi = opt.xi;
    
  for(int m=0; m<nidx; m++)
    {
      double force_x=0, force_y=0;
      xm[0] = x[idx[m]    ];
      xm[1] = x[idx[m]+N  ];
      for(int n=0; n<N; n++)
        {
	  rho2 = ( (xm[0]-x[n  ])*(xm[0]-x[n  ]) +
		   (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  if(rho2==0)
 	     continue;
	  force_x += q[n]*2.*(xm[0]-x[n]  )/rho2*(1-exp(-rho2*xi*xi));
	  force_y += q[n]*2.*(xm[1]-x[n+N])/rho2*(1-exp(-rho2*xi*xi));
        }
      force[idx[m]    ] = -force_x/opt.box[2];
      force[idx[m]+  N] = -force_y/opt.box[2];
      force[idx[m]+2*N] = 0;
    }  
}
#else
void SE1P_direct_k0(double* restrict phi,
        const int* restrict idx, int nidx,
        const double* restrict x,
        const double* restrict q, int N,
        const ewald_opts opt)
{
    double p;
    const double xi = opt.xi;
    double egamma = 0.57721566490153286061;

#ifdef _OPENMP
#pragma omp parallel for private(p)
#endif    
    for(int m=0; m<nidx; m++) {
      p=0;
      double xm[2] = {x[idx[m]    ],
		      x[idx[m]+N  ]};
      for(int n=0; n<N; n++) {
	double rho2 = ( (xm[0]-x[n  ])*(xm[0]-x[n  ]) +
			(xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	if(rho2>__DBL_EPSILON__)
	  p += -q[n]*(gsl_sf_expint_E1(rho2*xi*xi)+log(rho2*xi*xi)+egamma);
      }
      phi[m] = p/opt.box[2];
    }
}
#endif

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
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","K01P",N,num_eval);
    }

    // call kernel
    SE1P_direct_k0(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}
