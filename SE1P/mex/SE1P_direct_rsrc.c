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
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));

    opt->box[0] = box[0];

    // layers: mandatory for ewald sums that are truncated 
    const mxArray* mx_layers = mxGetField(mx_opt,0,"layers");
    if(mx_layers)
	opt->layers = (int)mxGetScalar(mx_layers);
    else
	opt->layers = -1;

    // rc: mandatory for short-range real sum 
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

#ifdef FORCE
void SE1P_direct_real_rc(double* restrict force, 
			 const int* restrict idx, int nidx,
			 const double* restrict x, 
			 const double* restrict q, int N,
			 const ewald_opts opt)
{
  double xi = opt.xi;
  double xi2 = xi*xi;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int m=0; m<nidx; m++)
    {
      double f[] = {0,0,0};
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p0 = -opt.layers; p0<=opt.layers; p0++)
	    {
	      if(idx[m] == n && p0 == 0)
		continue;

	      double rvp[] = {rvec[0]+p0*opt.box[0],rvec[1],rvec[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      double r2 = r*r;

	      if(r > opt.rc) continue;
	      
	      double c = qn*(2*xi/sqrt(PI)*exp(-xi2*r2)+ erfc(xi*r)/r)/r2;
	      f[0] += c*rvp[0];
	      f[1] += c*rvp[1];
	      f[2] += c*rvp[2];
	    }
	}
      force[m       ] = -f[0];
      force[m+  nidx] = -f[1];
      force[m+2*nidx] = -f[2];
    }
}
#else
void SE1P_direct_real_rc(double* restrict phi, 
			 const int* restrict idx, int nidx,
			 const double* restrict x, 
			 const double* restrict q, int N,
			 const ewald_opts opt)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m=0; m<nidx; m++)
    {
	double p = 0;
	for(int n=0; n<N; n++)
	{
	  double rvec[] = {x[idx[m]    ]-x[n    ],
			   x[idx[m]+N  ]-x[n+N  ],
			   x[idx[m]+2*N]-x[n+2*N]};
	    
	    double qn = q[n];

	    for(int p0 = -opt.layers; p0<=opt.layers; p0++)
	      {
		if(idx[m] == n && p0 == 0)
		  continue;
		
		double r = sqrt((rvec[0]+p0*opt.box[0])*
				(rvec[0]+p0*opt.box[0])+
				rvec[1]*rvec[1]+
				rvec[2]*rvec[2]
				);
		if(r < opt.rc)
		  p += qn*erfc(opt.xi*r)/r;
	      }
	}
	phi[m] += p;
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
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","RSRC1P",N,num_eval);
	mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
		  opt.xi,opt.rc,opt.layers);
    }

    // call kernel
    SE1P_direct_real_rc(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}

