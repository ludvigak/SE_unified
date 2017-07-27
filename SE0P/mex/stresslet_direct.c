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
      mexPrintf("[FS Stresslet Direct ] MEX N=%d ",N);

    // call kernel
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m=0; m<N; m++)
    {
	double p[] = {0,0,0};
	for(int n = 0; n<N; n++){
	  double  q[] = {f[n],     f[n+N],   f[n+2*N]};
	  double  s[] = {f[n+3*N], f[n+4*N], f[n+5*N]};
	  double  r[] = {x[m]-x[n], x[m+N]-x[n+N],x[m+2*N]-x[n+2*N]};
	  
	  double ri = 1.0/sqrt(dot(r,r));
	  double ri5 = ri*ri*ri*ri*ri;

	  double rdots = dot(r,s);
	  double rdotq = dot(r,q);
	  
	  double c = -6.0* rdots * rdotq * ri5;

	  if(m==n)
	      continue;
	  p[0] += c*r[0];
	  p[1] += c*r[1];
	  p[2] += c*r[2];
	}

	u[m    ] = p[0];
	u[m+  N] = p[1];
	u[m+2*N] = p[2];
      }
}
