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

static inline double dot(double x[], double y[])
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
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
  const double* restrict qvec = f;
  const double* restrict nvec = f + 3*M;

  
  // output
  U = mxCreateDoubleMatrix(M, 3, mxREAL);
  double* restrict u = mxGetPr(U);

  if(VERBOSE)
    mexPrintf("[FS Stresslet Real space ] MEX N=%d ",M);


  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  
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

    double p[3];
    p[0] = 0; p[1] = 0; p[2] = 0;

    // First element is the target itself; see MATLAB doc for rangesearch.
    for(int n = 1; n<N; n++){ 
      int idxn = (int) idx[n]-1;
      double   qn[] = {qvec[idxn], qvec[idxn+M], qvec[idxn+2*M]};
      double   nn[] = {nvec[idxn], nvec[idxn+M], nvec[idxn+2*M]};
      double rvec[] = {x[m]-x[idxn], x[m+M]-x[idxn+M],x[m+2*M]-x[idxn+2*M]};
      double      r = dis[n];
      double      r2= r*r;

      double c = xi2*r2;    
      double expc = exp(-c);
      // Hasimoto
      double C = -2/(r2*r2)*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(PI)*(3.0+2.0*c)*expc );
      double D = 4/sqrt(PI)*xi3*expc;

      double rdotn = dot(rvec, nn);
      double rdotq = dot(rvec, qn);
      double ndotq = dot(nn, qn);
      double Kr = C*rdotn*rdotq + D*ndotq;
      double Kn = D*rdotq;
      double Kq = D*rdotn;
      for (int i=0; i<3; i++)
      {
	  p[i] += Kr*rvec[i] + Kn*nn[i] + Kq*qn[i];
      }
    }
    u[m    ] = p[0];
    u[m+  M] = p[1];
    u[m+2*M] = p[2];
  }
}
