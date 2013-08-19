// Simple real space interaction computation that saves matrix for later use.
// Does not make use of interaction symmetries.

#include "mex_compat.h"
#include "math.h"

#define X    prhs[0] // Points (Nx3)
#define NN   prhs[1] // Point normals (Nx3)
#define IDX  prhs[2] // Target point indices
#define NBOX prhs[3] // Number of layers
#define RC   prhs[4] // Cutoff radius
#define XI   prhs[5] // Ewald parameter
#define BOX  prhs[6] // Ewald parameter

#define AMAT  plhs[0]  // Output (3x3 cell array of MxN matrices)

#define PI 3.141592653589793

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
    const double* restrict xvec = mxGetPr(X);
    const double* restrict nvec = mxGetPr(NN);
    const double* restrict box = mxGetPr(BOX);
    const int* restrict idx = mxGetData(IDX);

    // Get constants
    const int nbox = mxGetScalar(NBOX);
    const double rc = mxGetScalar(RC);
    const double xi = mxGetScalar(XI);

    // Get sizes
    const int nidx = mxGetNumberOfElements(IDX);

    const mwSize *Xdims;
    Xdims = mxGetDimensions(X);
    const int N = Xdims[0];

    // Output matrices
    int outDims[2];
    outDims[0] = nidx;
    outDims[1] = N;
    size_t malloc_tot = N*nidx*6*sizeof(double);
    ASSERT(malloc_tot < MALLOC_MAX, "MALLOC_MAX exceeded");
    AMAT = mxCreateCellMatrix(3,3);
    double* restrict amat[3][3];
    for(int i=0; i<3; i++)
	for(int j=i; j<3; j++)
	{
	    mxArray *matPt = mxCreateNumericArray(2,outDims,mxDOUBLE_CLASS,mxREAL);
	    amat[i][j] = mxGetPr( matPt );
	    mxSetCell(AMAT, i+3*j, matPt);
	}

    // Checking
    const mxClassID IDXclass =  mxGetClassID(IDX);
    if(IDXclass != mxINT32_CLASS)
	mexErrMsgTxt("Vector idx must be of type int32.");

#ifdef _OPENMP
#pragma omp parallel shared(amat,nvec,xvec,idx,box) default(none)
#endif
{
    
    double A[3][3];
    double AA[3][3];
    double xn[3];
    double nn[3];
    double r[3];
    double rm[3];
    double rsq,rcsq;
    int idxm,idxn;
    int m,n;
    int i1, i2, i3, k1, k2;

    rcsq = rc*rc;
    
    // for all points (sources)
#ifdef _OPENMP  
#pragma omp for
#endif
    for(n=0; n<N; n++)
    {
	idxn = n+1;
        // Put source vector in xs
        xn[0] = xvec[idxn    -1];
        xn[1] = xvec[idxn+N  -1];
        xn[2] = xvec[idxn+2*N-1];
        
        // Put normal vector in ns
        nn[0] = nvec[idxn    -1];
        nn[1] = nvec[idxn+N  -1];
        nn[2] = nvec[idxn+2*N-1];
        
        // for all evaluation points
        for(m=0; m<nidx; m++)               
        {
            // rm = xm - xn
	    idxm = idx[m];
            rm[0] = xvec[idxm    -1]-xn[0];  // indirect indexing OK in outer loop
            rm[1] = xvec[idxm+N  -1]-xn[1];
            rm[2] = xvec[idxm+2*N-1]-xn[2];
            
            // AA=0
            for(k1=0;k1<3;k1++)
                for(k2=0;k2<3;k2++)
                    AA[k1][k2]=0;
            
            // *** Collect contribution from all images in AA
            for(i1 = -nbox; i1<=nbox; i1++) { // image boxes
                for(i2 = -nbox; i2<=nbox; i2++) {
                    for(i3 = -nbox; i3<=nbox; i3++)
                    {
                        if(i1==0 && i2==0 && i3==0 && idxm==idxn) // skip self
                            continue;
                        r[0] = rm[0]+box[0]*i1;
                        r[1] = rm[1]+box[1]*i2;
                        r[2] = rm[2]+box[2]*i3;
                        rsq  = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
                        if(rsq <= rcsq)
                        {
                            op_A(A,r,nn,xi);
                            for(k1=0;k1<3;k1++)
                                for(k2=0;k2<3;k2++)
                                    AA[k1][k2]+=A[k1][k2];
                        }
                    }
                }
            }
            
	    // Set A{k1}{k2}[m][n] += AA[k1][k2]
            // {k1,k2}={0,1,2} over components of source and target, k1<=k2
            // A{}{} is nidx x N
	    // No risk of racing between threads here
            for(int k1=0; k1<=2; k1++)
		for(int k2=k1; k2<=2; k2++)
		    amat[k1][k2][m+nidx*n] += AA[k1][k2];
	}
    }
 }
}
