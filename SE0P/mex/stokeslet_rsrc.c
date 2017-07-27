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

inline double dot(double * a,double * b)
{
  return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
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

    double p[3], rdotf;
    p[0] = 0; p[1] = 0; p[2] = 0;

    // First element is the target itself; see MATLAB doc for rangesearch.
    for(int n = 1; n<N; n++){ 
      int idxn = (int) idx[n]-1;
      double   fn[] = {f[idxn],      f[idxn+M],        f[idxn+2*M]};
      double rvec[] = {x[m]-x[idxn], x[m+M]-x[idxn+M], x[m+2*M]-x[idxn+2*M]};
      double      r = dis[n];
      double      r2= r*r;
      double  xiexp = xi*exp(-xi*xi*r2);
      double     c1 = 2.0*( xiexp / (sqrt(PI)*r2) + erfc(xi*r) / (2*r*r2) );
      double     c2 = -4.0/sqrt(PI)*xiexp;

      rdotf = dot(rvec,fn);
      
      p[0] += (c1*r2 + c2)*fn[0] + c1*rdotf*rvec[0];
      p[1] += (c1*r2 + c2)*fn[1] + c1*rdotf*rvec[1];
      p[2] += (c1*r2 + c2)*fn[2] + c1*rdotf*rvec[2];
    }
    u[m    ] = p[0];
    u[m+  M] = p[1];
    u[m+2*M] = p[2];
  }
}
