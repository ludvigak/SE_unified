#include "mex.h"
#include "math.h"

#define P    prhs[0] // Evaluation points (Matlab/Fortran 1-based idices)
#define X    prhs[1] // Source locations
#define F    prhs[2] // Source strengths
#define XI   prhs[3] // Ewald parameters

#define PHI plhs[0]  // Output

#define PI 3.141592653589793

// select real space scaling tensor
#ifdef HASIMOTO
#include "hasimoto_op_real.h"
#elif BEENAKKER
#include "beenakker_op_real.h"
#else
#error "Must provide -D<method> to compiler"
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    if(VERBOSE)
	mexPrintf("Stokes Ewald kernel (real): %s\n", OP_TAG);

    // input dims
    const int N = mxGetN(X);
    if ( mxGetM(F) != 3 || mxGetM(X) != 3 || mxGetN(F) != N)
    {
	mexPrintf("Stokes Ewald kernel: ERROR\n\t X and F must be 3-by-N\n");
	return;
    }

    const double xi = mxGetScalar(XI);
    
    const double* restrict p   = mxGetPr(P);
    const double* restrict x   = mxGetPr(X);
    const double* restrict f   = mxGetPr(F);

    PHI = mxCreateDoubleMatrix(3, 1, mxREAL);
    double* phi = mxGetPr(PHI);

    double A[3][3];
    double r[3];
    int j;
    for(j=0; j<3*N; j+=3)            // for all particles
    {
	r[0] = p[0] - x[j  ];
	r[1] = p[1] - x[j+1];
	r[2] = p[2] - x[j+2];

	op_A(A,r,xi);                              // phi += A*f
	phi[0] += A[0][0]*f[j] + A[0][1]*f[j+1] + A[0][2]*f[j+2];
	phi[1] += A[1][0]*f[j] + A[1][1]*f[j+1] + A[1][2]*f[j+2];
	phi[2] += A[2][0]*f[j] + A[2][1]*f[j+1] + A[2][2]*f[j+2];
    }
}
