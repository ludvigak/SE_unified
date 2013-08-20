#include "mex.h"
#include <math.h>

#define X  prhs[0]
#define Q  prhs[1]
#define XI prhs[2]
#define Z  prhs[3]

#define H plhs[0]  // Output

#define SQRT_PI 1.772453850905516

/* no input checking is done */
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    // input dims
    const int N = mxGetM(X);
    const int M = mxGetM(Z); 
    const double* x = mxGetPr(X);
    const double* q = mxGetPr(Q);
    const double* z = mxGetPr(Z);

    H = mxCreateDoubleMatrix(M, 1, mxREAL);
    double* restrict h = mxGetPr(H);

    double xi = mxGetScalar(XI);
    double xi2 = xi*xi;

    double z_ij,p;
    for(int i=0; i<M; i++)
    {
	p=0;
	for(int j=0; j<N; j++)
	{
	    z_ij = x[j+2*N]-z[i];
	    p += q[j]*(exp(-xi2*z_ij*z_ij)/xi + SQRT_PI*z_ij*erf(xi*z_ij));
	}
	h[i] = p;
    }    
}
