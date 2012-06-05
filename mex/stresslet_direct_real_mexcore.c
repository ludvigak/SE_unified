#include "mex.h"
#include "math.h"

#define X    prhs[0] // Target points (Nx3)
#define IDX  prhs[1] // Target point indices
#define XN   prhs[2] // Source point (1x3)
#define NN   prhs[3] // Source point normal (1x3)
#define IDXN prhs[4] // Source point index
#define NBOX prhs[5] // Number of layers
#define XI   prhs[6] // Ewald parameter
#define BOX  prhs[7] // Ewald parameter

#define TMP  plhs[0]  // Output (Mx3)

#define PI 3.141592653589793

// select k-space scaling tensor
#ifdef BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    if(VERBOSE)
		mexPrintf("Stresslet real space core: %s\n", OP_TAG);

	// Get input arrays
    const double* restrict x = mxGetPr(X);
    const double* restrict box = mxGetPr(BOX);
	const int* restrict idx = mxGetData(IDX);
    const double* restrict xn = mxGetPr(XN);
    const double* restrict nn = mxGetPr(NN);

	// Get constants
    const int idxn = mxGetScalar(IDXN);
    const int nbox = mxGetScalar(NBOX);
    const double xi = mxGetScalar(XI);

	// Get sizes
	const int nidx = mxGetNumberOfElements(IDX);

    const mwSize *Xdims;
    Xdims = mxGetDimensions(X);
	const int N = Xdims[0];


    // output array
	int tmpDims[2];
	tmpDims[0] = 3*nidx;
	tmpDims[1] = 3;

	// Checking
	const mxClassID IDXclass =  mxGetClassID(IDX);
	if(IDXclass != mxINT32_CLASS)
		mexErrMsgTxt("Vector idx must be of type int32.");

	TMP = mxCreateNumericArray(2,tmpDims,mxDOUBLE_CLASS,mxREAL);
    double* restrict tmp = mxGetPr( TMP );

    double A[3][3];
	double n[3];
	double r[3];
	double rm[3];
    int i1, i2, i3, m, k1, k2;

	n[0] = nn[0];
	n[1] = nn[1];
	n[2] = nn[2];

    for(m=0; m<nidx; m++)                          // for all evaluation points
    {

	// rm = xm - xn
	rm[0] = x[idx[m]    -1]-xn[0];             // indirect indexing OK in outer loop
	rm[1] = x[idx[m]+N  -1]-xn[1];
	rm[2] = x[idx[m]+2*N-1]-xn[2];

	for(i1 = -nbox; i1<=nbox; i1++)                           // image boxes
	    for(i2 = -nbox; i2<=nbox; i2++)
			for(i3 = -nbox; i3<=nbox; i3++)
			{
				if(i1==0 && i2==0 && i3==0 && idx[m]==idxn)    // skip self
			    	continue;

				r[0] = rm[0]+box[0]*i1;
				r[1] = rm[1]+box[1]*i2;
				r[2] = rm[2]+box[2]*i3;

				op_A(A,r,n,xi); // Get A, but only use lower part

				// Set tmp[m+k1*N][k2] = tmp[m+k1*N][k2] + A[k1][k2],
				// {k1,k2}={0,1,2} over components of source and target
				// tmp is 3Nx3, so tmp[m+k1*N][k2]=tmp[m+N(k1+3k2)]
				tmp[m+nidx*0] = tmp[m+nidx*0] + A[0][0];
				tmp[m+nidx*1] = tmp[m+nidx*1] + A[0][1];
				tmp[m+nidx*2] = tmp[m+nidx*2] + A[0][2];
				tmp[m+nidx*3] = tmp[m+nidx*3] + A[0][1];
				tmp[m+nidx*4] = tmp[m+nidx*4] + A[1][1];
				tmp[m+nidx*5] = tmp[m+nidx*5] + A[1][2];
				tmp[m+nidx*6] = tmp[m+nidx*6] + A[0][2];
				tmp[m+nidx*7] = tmp[m+nidx*7] + A[1][2];
				tmp[m+nidx*8] = tmp[m+nidx*8] + A[2][2];
			}
    }
}
