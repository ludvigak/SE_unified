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

    // layers: mandatory for ewald sums that are truncated 
    const mxArray* mx_layers = mxGetField(mx_opt,0,"layers");
    if(mx_layers)
	opt->layers = (int)mxGetScalar(mx_layers);
    else
	opt->layers = -1;
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
}

#ifdef FORCE
void SE1P_direct_fd(double* restrict force, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
  double f[3];
  const double xi   = opt.xi;
  double xi2        = xi*xi;
  double TwoPiOverL = 2.*PI/opt.box[2];
  /* get ready for gsl integration*/
   gsl_integration_workspace * w 
     = gsl_integration_workspace_alloc (STACK_SIZE);

  for(int m=0; m<nidx; m++)
    {
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};

      f[0] = 0; f[1] = 0; f[2] = 0;
      for(int n = 0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n],xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double rho2 = rvec[0]*rvec[0] + rvec[1]*rvec[1];
	  double b   = rho2*xi2;
	  double qn  = q[n];
	  for(int j2 = -opt.layers; j2<=opt.layers; j2++)
	    {
	      if(j2 == 0)
		continue;
	      
	      double k3  = TwoPiOverL*j2;
	      double k3z = -k3*rvec[2];
	      
	      double a   = k3*k3/(4.*xi2);
	      //K0  = IncompBesselK0_Simpson(1e-15, &rep, a,b,1);
	      double K0  = call_gsl_bessel_integrator(a,b,w,1);
	      f[0] += 2.*qn*xi2*cos(k3z)*rvec[0]*K0;
	      f[1] += 2.*qn*xi2*cos(k3z)*rvec[1]*K0;
	      //K0  = IncompBesselK0_Simpson(1e-15, &rep, a,b,0);
	      K0  = call_gsl_bessel_integrator(a,b,w,0);
	      f[2] += -qn*k3*sin(k3z)*K0;
	    }
	}
      force[idx[m]    ] = -f[0]/(opt.box[2]);
      force[idx[m]+  N] = -f[1]/(opt.box[2]);
      force[idx[m]+2*N] = -f[2]/(opt.box[2]);
    }

    gsl_integration_workspace_free (w);
}
#else
void SE1P_direct_fd(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
  
  double p;
  const double xi   = opt.xi;
  double xi2        = xi*xi;
  double TwoPiOverL = 2.*PI/opt.box[2];
  //  int rep;

#ifdef _OPENMP
#pragma omp parallel for private(p)
#endif
  for(int m=0; m<nidx; m++)
    {
      double xm[3] = {x[idx[m]    ],
		      x[idx[m]+N  ],
		      x[idx[m]+2*N]};
      p = 0;
      for(int j2 = 1; j2<=opt.layers; j2++) {
	double k3  = TwoPiOverL*j2;
	
	double a   = k3*k3/(4.*xi2);
	for(int n = 0; n<N; n++) {
	  double z   = xm[2]-x[n+2*N];
	  double rho2= ( (xm[0]-x[n  ])*(xm[0]-x[n  ])+
			 (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  double b   = rho2*xi2;
	  double qn  = q[n];

	  //	      K0  = IncompBesselK0_Simpson(1e-15 ,&rep, a, b,0);
	  double K0 = computeK0(a,b);
	  double k3z = -k3*z;
	  p += 2*qn*cos(k3z)*K0;
	  
	  
	}
      }
      phi[m] = p/(opt.box[2]);
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
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","FD1P",N,num_eval);
	mexPrintf("xi = %.2f, layers=%d\n", opt.xi,opt.layers);
    }

    // call kernel
    SE1P_direct_fd(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}
