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
void SE1P_direct(double* restrict force, 
		 const int* restrict idx, int nidx,
		 const double* restrict x, 
		 const double* restrict q, int N,
		 const ewald_opts opt)
{
  double f[3];
#ifdef _OPENMP
#pragma omp parallel for private(f)
#endif
  for(int m=0; m<nidx; m++)
    {
      f[0] = 0; f[1] = 0; f[2] = 0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p2 = -opt.layers; p2<=opt.layers; p2++)
	    {
	      if(idx[m] == n && p2 == 0)
		continue;

	      double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      double r3 = r*r*r;
	      
	      double c = qn/r3;
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
void SE1P_direct(double* restrict phi, 
		 const int* restrict idx, int nidx,
		 const double* restrict x, 
		 const double* restrict q, int N,
		 const ewald_opts opt)
{
  double p;
#ifdef _OPENMP
#pragma omp parallel for private(p)
#endif
  for(int m=0; m<nidx; m++)
    {
      p=0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int p2 = -opt.layers; p2<=opt.layers; p2++) {
	for(int n=0; n<N; n++)
	  {
	    double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	    double qn = q[n];
	    
	    if(idx[m] == n && p2 == 0)
	      continue;
	    
	    double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	    double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	    
	    p += qn/r;
	  }
      }
      phi[m] = p;
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
      mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","D1P",N,num_eval);

    // call kernel
    SE1P_direct(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}
