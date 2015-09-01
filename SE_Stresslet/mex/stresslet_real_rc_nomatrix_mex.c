#include "stresslet_real_rc.h"

#define X    prhs[0] // Points (Nx3)
#define NVEC prhs[1] //        (Nx3)
#define FVEC prhs[2] //        (Nx3)
#define BOX  prhs[3] // Box size
#define RC   prhs[4] // Cutoff radius
#define XI   prhs[5] // Ewald parameter

#define PHI plhs[0]

//==============================================
// ==== MAIN MEX FUNCTION
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Setup variables
    double* restrict phi; 
   
    // Get input arrays
    const double* restrict x = mxGetPr(X);
    const double* restrict nvec = mxGetPr(NVEC);
    const double* restrict fvec = mxGetPr(FVEC);
    const double* restrict box = mxGetPr(BOX);

    // Get constants
    const double xi = mxGetScalar(XI);
    const double rc = mxGetScalar(RC);

    // Get number of particles
    const int N = mxGetM(X);

    // Prepare output
    PHI = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);

    // Compute
    compute_rsrc_direct (x, nvec, fvec, N, box, xi, rc, &phi);

    // Return
    mxSetM(PHI, N);
    mxSetN(PHI, 3);
    mxSetPr(PHI, phi);
}


