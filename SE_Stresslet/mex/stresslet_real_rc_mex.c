#include "stresslet_real_rc.h"

#define X    prhs[0] // Points (Nx3)
#define NVEC prhs[1] // Points (Nx3)
#define BOX  prhs[2] // Box size
#define RC   prhs[3] // Cutoff radius
#define XI   prhs[4] // Ewald parameter

#define AMAT plhs[0]
#define ROW plhs[1]
#define COL plhs[2]
#define VAL plhs[3]
#define PER plhs[4]

//==============================================
// ==== MAIN MEX FUNCTION
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Setup variables
    int i,j,numel;

    int* restrict row;
    int* restrict col;
    int* restrict idx_in_array;
    int* restrict buck_size;
    double* restrict val[3][3];

    struct timeval tic, toc;
    double time_spent;
   
    // Get input arrays
    const double* restrict x = mxGetPr(X);
    const double* restrict nvec = mxGetPr(NVEC);
    const double* restrict box = mxGetPr(BOX);

    // Get constants
    const double xi = mxGetScalar(XI);
    const double rc = mxGetScalar(RC);

    // Get number of particles
    const int N = mxGetM(X);

    //===============================================
    // BUILD TRIPLETS
    //
    gettimeofday(&tic, NULL);

    get_rs_triplets ( x, nvec, N, 
		      box, xi, rc, nlhs,
		      &row, &col, val, &buck_size, &idx_in_array, &numel
		      );

    gettimeofday(&toc,NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
	mexPrintf("[RSRC] = Total triplet time: %.3f seconds.\n", time_spent);

    //===============================================
    // BUILD SPARSE MATRICES
    //
    gettimeofday(&tic, NULL);
    AMAT = mxCreateCellMatrix(3,3);
    mxArray *matPt;
    mwIndex* restrict jc[6];
    mwIndex* restrict ir[6];
    int k1,k2;

    struct timeval tic2, toc2;

    gettimeofday(&tic2, NULL);
    // Setup row and column lists
    for(j=0; j<6; j++)
    {
	jc[j] = mxMalloc((N+1)*sizeof(mwIndex));
	ir[j] = mxMalloc(numel*sizeof(mwIndex));
    }
    jc[0][0] = 0;
    // Setup first column list with bucket sort results
    for(j=0; j<N; j++)
    {
	jc[0][j+1] = jc[0][j]+buck_size[j];
    }
    // Copy to remaining columns
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(j=0; j<N+1; j++)
    {
	i = jc[0][j];
	jc[1][j] = i;
	jc[2][j] = i;
	jc[3][j] = i;
	jc[4][j] = i;
	jc[5][j] = i;
    }
    // Copy row to all rows
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i=0; i<numel; i++)
    {
	j = row[i];
	ir[0][i] = j;
	ir[1][i] = j;
	ir[2][i] = j;
	ir[3][i] = j;
	ir[4][i] = j;
	ir[5][i] = j;
    }
    gettimeofday(&toc2, NULL);
    if(VERBOSE)
	mexPrintf("[RSRC] Rows/cols copied in %.3f seconds.\n", DELTA(tic2,toc2));



    mxFree(buck_size); // Free final bucket list
    if(nlhs==1)
	mxFree(row); // Free row list if only returning matrix

    int k1_l[] = {0,     0,     0,     1,     1,     2};
    int k2_l[] = {0,     1,     2,     1,     2,     2};

    double* restrict pr[6];

    // Setup vectors for matrix data and permute to sorted
    gettimeofday(&tic2, NULL);
    for(j=0; j<6; j++)
    {		
	k1 = k1_l[j];
	k2 = k2_l[j];
	pr[j]  = mxMalloc(numel*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
	for(i=0;i<numel;i++)
	{	
	    pr[j][i] = val[k1][k2][ idx_in_array[i] ];
	}

	// Free original vals if no longer needed
	if(nlhs==1)
	    mxFree(val[k1][k2]);
    }
    gettimeofday(&toc2, NULL);
    if(VERBOSE)
	mexPrintf("[RSRC] Data sorted in %.3f seconds.\n", DELTA(tic2,toc2));

    // Add sparse matrices to cell array
    for(j=0; j<6; j++)
    {	
	k1 = k1_l[j];
	k2 = k2_l[j];
	matPt = mxCreateSparse(0,0,0,mxREAL);
	mxSetPr(matPt,pr[j]);
	mxSetJc(matPt,jc[j]);
	mxSetIr(matPt,ir[j]);

	mxSetM(matPt,N);
	mxSetN(matPt,N);
	mxSetNzmax(matPt,numel);
	mxSetCell(AMAT, k1 + 3*k2, matPt);
    }

    gettimeofday(&toc,NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
	mexPrintf("[RSRC] = Matrices assembled in %.3f seconds.\n", time_spent);

    //====================================================
    // (optional) OUTPUT TRIPLETS AND PERMUTATION LIST
    //
    if(nlhs>1)
    {
	gettimeofday(&tic,NULL);

	// Convert to MATLAB indexing
	for(i=0;i<numel;i++)
	    col[i]++;
	for(i=0;i<numel;i++)
	    row[i]++;
	for(i=0;i<numel;i++)
	    idx_in_array[i]++;

	COL = mxCreateNumericArray(0,0,mxUINT32_CLASS,mxREAL);
	ROW = mxCreateNumericArray(0,0,mxUINT32_CLASS,mxREAL);
	VAL = mxCreateCellMatrix(3,3);
   	
	mxArray *arrPt;
	for(i=0;i<=2;i++)
	    for(j=i;j<=2;j++)
	    {
		arrPt = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
		mxSetPr(arrPt, val[i][j]);
		mxSetM(arrPt, numel);
		mxSetN(arrPt, 1);
		mxSetCell(VAL, i + 3*j, arrPt);
	    }

	mxSetData(COL, col);
	mxSetData(ROW, row);

	mxSetM(COL, numel);
	mxSetM(ROW, numel);

	mxSetN(COL, 1);
	mxSetN(ROW, 1);


	PER = mxCreateNumericArray(0,0,mxINT32_CLASS,mxREAL);
	mxSetData(PER, idx_in_array);
	mxSetN(PER,1);
	mxSetM(PER,numel);

	gettimeofday(&toc,NULL);
	time_spent = DELTA(tic,toc);
	if(VERBOSE)
	    mexPrintf("[RSRC] Extra output prepared in %.3f seconds.\n", time_spent);
    }
}


