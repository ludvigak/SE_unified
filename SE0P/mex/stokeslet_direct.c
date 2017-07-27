#include "mex.h"
#include "math.h"

#define X   prhs[0] // Source locations
#define F   prhs[1] // Source strengths

#define U   plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define PI 3.141592653589793

inline double dot(double * a,double * b)
{
  return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

/* no input checking is done */
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    // input dims
    const int N = mxGetM(X);
   
    const double* restrict x = mxGetPr(X);
    const double* restrict f = mxGetPr(F);

    U = mxCreateDoubleMatrix(N, 3, mxREAL);
    double* restrict u = mxGetPr(U);

    if(VERBOSE)
      mexPrintf("[FS Stokeslet Direct ] MEX N=%d ",N);

    // call kernel
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m=0; m<N; m++)
    {
	double p[] = {0,0,0};
	for(int n = 0; n<N; n++){
	  double fn[] = {f[n], f[n+N], f[n+2*N]};
	  double  r[] = {x[m]-x[n], x[m+N]-x[n+N],x[m+2*N]-x[n+2*N]};
	  
	  double rdotf = dot(r,fn);
	  double ri = 1.0/sqrt(dot(r,r));
	  double ri3 = ri*ri*ri;

	  if(m==n)
	    continue;
	  p[0] += ri*fn[0] + rdotf*ri3*r[0];
	  p[1] += ri*fn[1] + rdotf*ri3*r[1];
	  p[2] += ri*fn[2] + rdotf*ri3*r[2];
	}
	u[m    ] = p[0];
	u[m+  N] = p[1];
	u[m+2*N] = p[2];
      }
}
