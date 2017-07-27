#include "mex.h"
#include "math.h"

#define X   prhs[0] // Source locations
#define F   prhs[1] // Source strengths
#define IDX prhs[2] // Source indeces
#define DIS prhs[3] // Distances
#define XI  prhs[4] // Ewald Param

#define U   plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define PI 3.141592653589793

inline void cross(double * a,double * b, double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{

  // input target
  const int M = mxGetM(X);
  const double xi = (double) mxGetScalar(XI);
  const double* restrict x = mxGetPr(X);
  const double* restrict f = mxGetPr(F);
  
  // output
  U = mxCreateDoubleMatrix(M, 3, mxREAL);
  double* restrict u = mxGetPr(U);

  if(VERBOSE)
    mexPrintf("[FS Rotlet Real space ] MEX N=%d ",M);  
  
  // Loop through the cell 
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int m=0; m<M; m++)  {
    const mxArray * _IDX = mxGetCell(IDX, m);
    const mxArray * _DIS = mxGetCell(DIS, m);
    const int N= mxGetN(_IDX);  // number of the source points in nblist
    double* idx= (double*) mxGetPr(_IDX);// index pointer of the source points in nblist
    double* dis= mxGetPr(_DIS); // distance pointer of the source points in nblist

    double p[3], fxr[3];
    p[0] = 0; p[1] = 0; p[2] = 0;

    // First element is the target itself; see MATLAB doc for rangesearch.
    for(int n = 1; n<N; n++){ 
      int idxn = (int) idx[n]-1;
      double   fn[] = {f[idxn], f[idxn+M], f[idxn+2*M]};
      double rvec[] = {x[m]-x[idxn], x[m+M]-x[idxn+M],x[m+2*M]-x[idxn+2*M]};
      double      r = dis[n];
      double      r2= r*r;

      cross(fn,rvec,fxr);
      double A = ( erfc(xi*r)/r + 2.0*xi*exp(-xi*xi*r2)/sqrt(PI) )/r2;

      p[0] += A*fxr[0];
      p[1] += A*fxr[1];
      p[2] += A*fxr[2];
    }
    u[m    ] = p[0];
    u[m+  M] = p[1];
    u[m+2*M] = p[2];
  }
}
