#include "mex.h"
#include "math.h"
#include "string.h"
#include "time.h"

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


#define PI 3.141592653589793

#ifndef VERBOSE
#define VERBOSE 0
#endif

#ifdef BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

// ==== STRUCT AND COMPARISON FUNCTION FOR SORTING
typedef struct
{
    int idx_in_array, row; 
} idx_qsort_buck;

static int compare_qsort_buck(const void* a, const void* b)
{
    int temp = ( * ((idx_qsort_buck*) a) ).row - 
	( * ((idx_qsort_buck*) b) ).row;
    if (temp > 0)
	return 1;
    else if (temp < 0)
	return -1;
    else
	return 0;
}

// ==== MAIN MEX FUNCTION
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Setup variables
    int i,j,k,idx_s,idx_t,ip;
    int head_idx, N, ncell_tot;
    int ncell[3], icell[3], home_cell[3];
    int* restrict ll;
    int* restrict head;
    double boxmin, rsq, rcsq;
    double pshift[3], xs[3], xt[3], ns[3], nt[3], xr[3];
    double A1[3][3], A2[3][3];

    clock_t begin, end;
    double time_spent;
    begin = clock();

    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
   
    // Get input arrays
    const double* restrict x = mxGetPr(X);
    const double* restrict nvec = mxGetPr(NVEC);
    const double* restrict box = mxGetPr(BOX);

    // Get constants
    const double xi = mxGetScalar(XI);
    const double rc = mxGetScalar(RC);
    rcsq = rc*rc;

    // Get number of particles
    N = mxGetM(X);

    //==============================================================
    // BUILD CELL LIST
    //
    // TODO: Add some assertions to make sure rc not too big,
    // and that box can be divided into square cells.

    // Setup cell partitioning
    boxmin = box[0];
    if(box[1]<boxmin)
	boxmin = box[1];
    if (box[2]<boxmin)
	boxmin = box[2];
    const double rn = boxmin / floor(boxmin/rc);
    for(i=0;i<3;i++)
	ncell[i] = round( box[i]/rn );
    ncell_tot = ncell[0]*ncell[1]*ncell[2];

    if(VERBOSE)
    {
	mexPrintf("[RSRC] %s, xi=%g\n", OP_TAG, xi);
	mexPrintf("[RSRC] rc=%.3f, rn=%.3f\n", rc, rn);
	mexPrintf("[RSRC] box=(%g,%g,%g), ncell=(%d,%d,%d)\n", 
		  box[0],box[1],box[2],
		  ncell[0],ncell[1],ncell[2]);
    }

    // Prepare cell list
    ll = mxCalloc(N, sizeof(int));
    head = mxCalloc(ncell_tot, sizeof(int));
    for(i=0; i<ncell_tot; i++)
    	head[i] = -1;
    // Do cell partitioning 
    for(i=0; i<N; i++)
    {
	for(j=0; j<3; j++)
	    icell[j] = floor( x[i+N*j]/rn );
	head_idx = 
	    icell[0] +
	    icell[1]*ncell[0] + 
	    icell[2]*ncell[1]*ncell[0];
	ll[i] = head[head_idx];
	head[head_idx] = i;
    }

    //============================================================
    // CALCULATE INTERACTIONS
    //
    // For all vectors, go through neighbors and save interactions
    // in vectors that are used to create a sparse matrix
    // Allocate a guess based on average density +50%
    int maxel = round(  1.5 * N*N*4*PI*rc*rc*rc/3/(box[0]*box[1]*box[2]) );
    int numel = 0;
    int* restrict row = mxMalloc(maxel*sizeof(int));
    int* restrict col = mxMalloc(maxel*sizeof(int));
    double* restrict val[3][3];
    for(i=0;i<=2;i++)
	for(j=i;j<=2;j++)
	{
	    val[i][j] = mxMalloc(maxel*sizeof(double));
	}
    int compnum = 0; // Count distance comparisons
    for(idx_s=0;idx_s<N;idx_s++)
    {
	// Source point
	xs[0] = x[idx_s    ];
	xs[1] = x[idx_s+N  ];
	xs[2] = x[idx_s+2*N];
	// Source point normal vector
	ns[0] = nvec[idx_s    ];
	ns[1] = nvec[idx_s+N  ];
	ns[2] = nvec[idx_s+2*N];

	// Determine home cell
	for(j=0; j<3; j++)
	    home_cell[j] = floor( xs[j]/rn );	
	// Iterate through near cells (including home cell)
	for(ip=0; ip<27; ip++) 
	{
	    // Get neigh cell
	    icell[0] = home_cell[0] + px[ip];
	    icell[1] = home_cell[1] + py[ip];
	    icell[2] = home_cell[2] + pz[ip];
	    // Periodic wrap
	    for(j=0; j<3; j++)
	    {
		// (Could do this with mod)
		pshift[j] = 0;
		if(icell[j] >= ncell[j])
		{
		    icell[j] = 0;
		    pshift[j] = box[j];
		}
		else if(icell[j]<0)
		{
		    icell[j] = ncell[j]-1;
		    pshift[j] = -box[j];
		}
	    }
	    head_idx = 
		icell[0] +
		icell[1]*ncell[0] + 
		icell[2]*ncell[1]*ncell[0];	    
	    // Go through cell list
	    idx_t = head[head_idx];
	    while( idx_t != -1)
	    {
		if(idx_t > idx_s)
		{
		    // r points from s to t
		    for(j=0; j<3; j++)
			xr[j] = x[idx_t+j*N] + pshift[j] - xs[j];

		    // Check if we are within truncation radius
		    rsq = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
		    if(rsq <= rcsq)
		    {		
			// Source point normal vector
			nt[0] = nvec[idx_t    ];
			nt[1] = nvec[idx_t+N  ];
			nt[2] = nvec[idx_t+2*N];

			op_A_symm(A1,A2,xr,ns,nt,xi);

			/*
			op_A(A1,xr,ns,xi);
			xr[0] = -xr[0];
			xr[1] = -xr[1];
			xr[2] = -xr[2];
			op_A(A2,xr,nt,xi);
			*/

			compnum+=2;
			// Check if we still have allocated space
			if(numel>=maxel-1)
			{
			    mexPrintf("Reallocating...\n");
			    // Else, allocate 20% more
			    maxel = ceil(1.2*maxel);

			    row = mxRealloc(row, maxel*sizeof(int));
			    col = mxRealloc(col, maxel*sizeof(int));
			    for(i=0;i<=2;i++)
				for(j=i;j<=2;j++)
				    val[i][j] = mxRealloc(val[i][j], maxel*sizeof(double));
			}

			row[numel] = idx_t;
			col[numel] = idx_s;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][numel] = A1[i][j];
			    }
			numel++;

			row[numel] = idx_s;
			col[numel] = idx_t;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][numel] = A2[i][j];
			    }
			numel++;
		    }
		    
		}
		// Goto next point in cell chain
		idx_t = ll[idx_t];
	    } // End of neighbours in this cell
	} // End of cells
    } // End of particles

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(VERBOSE)
    {
	mexPrintf("[RSRC] Triplets generated in %.3f seconds.\n", time_spent);
    }

    //============================================
    // SORT RESULTS WITH BUCKET + QUICK SORT
    // Bucket sort on columns, then quicksort on rows
    // in each column
    begin = clock();
    int* restrict buck_size = mxMalloc(N*sizeof(int));
    int* restrict buck_count = mxMalloc(N*sizeof(int));
    int* restrict buck_pos = mxMalloc(N*sizeof(int));
    int* restrict idx_in_array = mxMalloc(numel*sizeof(int));
    int buck_idx, buck_max;

    // Init lists
    for(i=0;i<N;i++)
    {
	buck_size[i]=0;
	buck_count[i]=0;
    }

    // Count number of elements in each bucket (column)
    for(i=0;i<numel;i++)
    {
	buck_idx = col[i];
	buck_size[buck_idx]++;
    }

    // Cumulative addition to get locations of each bucket after sort,
    // + save largest bucket size for later.
    buck_pos[0] = 0;
    buck_max = 0;
    for(i=1;i<N;i++)
    {
	buck_pos[i] = buck_pos[i-1]+buck_size[i-1];
	if(buck_size[i-1]>buck_max)
	    buck_max = buck_size[i-1];
    }

    // Assign each element to a bucket, store permutations in idx_in_array
    for(i=0;i<numel;i++)
    { 
	buck_idx = col[i];
	idx_in_array[ buck_pos[buck_idx] + buck_count[buck_idx] ] = i;
	buck_count[buck_idx]++;
    }

    mxFree(buck_count); // Free counter
    if(nlhs==1)
	mxFree(col); // Free column list if only returning matrix
 
    // Quicksort on buckets (rows indices for a column)
    // Allocate list of pairs to fit largest bucket
    idx_qsort_buck* restrict qlist = mxMalloc(buck_max*sizeof(idx_qsort_buck));
    for(buck_idx=0;buck_idx<N;buck_idx++)
    {
	int begin = buck_pos[buck_idx];
	int size = buck_size[buck_idx];	    
	// Put (idx,row) pairs in list to be qsorted.
	for(i=0; i<size; i++)
	{
	    qlist[i].idx_in_array = idx_in_array[begin+i];
	    qlist[i].row = row[ idx_in_array[begin+i] ];
	}
	qsort(qlist, size, sizeof(idx_qsort_buck), compare_qsort_buck);	    
	for(i=0; i<size; i++)
	    idx_in_array[begin+i] = qlist[i].idx_in_array;
    }       

    mxFree(qlist); // Free temp list
    mxFree(buck_pos); // Free bucket list

    /* Removed insertion sort, dont think it beats qsort even on small buckets
    // Insertion sort on buckets
    double rowtmp;
    int idxtmp;
    for(buck_idx=0;buck_idx<N;buck_idx++)
    {
	int begin = buck_pos[buck_idx];
	int size = buck_size[buck_idx];

	for(i=begin+1;i<begin+size;i++)
	{
	    for(j=i; j>begin && row[ idx_in_array[j] ]<row[ idx_in_array[j-1] ]; j--)
	    {		    
		idxtmp = idx_in_array[j];
		idx_in_array[j] = idx_in_array[j-1];
		idx_in_array[j-1] = idxtmp;
	    }
	}
    }
    */


    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(VERBOSE)
	mexPrintf("[RSRC] Triplets sorted in %.3f seconds.\n", time_spent);


    //===============================================
    // BUILD SPARSE MATRICES
    //
    begin = clock();
    AMAT = mxCreateCellMatrix(3,3);
    mxArray *matPt;
    int k1,k2;
    mwIndex* restrict jc_t;
    mwIndex* restrict ir_t;
    mwIndex* restrict jc = mxMalloc((N+1)*sizeof(mwIndex));
    mwIndex* restrict ir = mxMalloc(numel*sizeof(mwIndex));

    jc[0] = 0;
    for(j=0; j<N; j++)
    {
	// Setup columns list with bucket sort results
	jc[j+1] = jc[j]+buck_size[j];
    }
    for(i=0; i<numel; i++)
    {
	// Reorder rows with permutation list
	ir[i] = row[ idx_in_array[i] ];
    }

    mxFree(buck_size); // Free final bucket list
    if(nlhs==1)
	mxFree(row); // Free row list if only returning matrix

    for(k1=0;k1<=2;k1++)
	for(k2=k1;k2<=2;k2++)
	{		
	    matPt = mxCreateSparse(0,0,0,mxREAL);
	    double* restrict pr = mxMalloc(numel*sizeof(double));
	    for(i=0; i<numel; i++)
	    {		
		// Reorder values with permutation list
		pr[i] = val[k1][k2][ idx_in_array[i] ];
	    }

	    if(nlhs==1)
		mxFree(val[k1][k2]); // Values no longer needed if only returning matrix

	    mxSetPr(matPt,pr);
	    if(k1==0 && k2==0)
	    {
		mxSetJc(matPt,jc);
		mxSetIr(matPt,ir);
	    }
	    else
	    {
		jc_t = mxMalloc((N+1)*sizeof(mwIndex));
		ir_t = mxMalloc(numel*sizeof(mwIndex));
		memcpy(jc_t,jc,(N+1)*sizeof(mwIndex));
		memcpy(ir_t,ir,numel*sizeof(mwIndex));
		mxSetJc(matPt,jc_t);
		mxSetIr(matPt,ir_t);
	    }
	    mxSetM(matPt,N);
	    mxSetN(matPt,N);
	    mxSetNzmax(matPt,numel);
	    mxSetCell(AMAT, k1 + 3*k2, matPt);
	}
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(VERBOSE)
	mexPrintf("[RSRC] Matrices assembled in %.3f seconds.\n", time_spent);
      
    //====================================================
    // (optional) OUTPUT TRIPLETS AND PERMUTATION LIST
    //
    if(nlhs>1)
    {
	begin = clock();

	for(i=0;i<numel;i++)
	{
	    // Convert to MATLAB indexing
	    col[i]++;
	    row[i]++;
	    idx_in_array[i]++;
	}

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

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	if(VERBOSE)
	    mexPrintf("[RSRC] Extra output prepared in %.3f seconds.\n", time_spent);
    }
}


