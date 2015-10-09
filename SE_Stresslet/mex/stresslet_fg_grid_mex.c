#include "mex.h"
#include "SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define XVEC   prhs[0] 
#define NVEC   prhs[1] 
#define FVEC   prhs[2] 
#define OPT prhs[3]

#define H_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

// Number of threads to use for the 9 units of work.
// Three threads has turned to be a sweet spot (locally).
#define NUM_THREADS 3

// Use the SSE code available
#define SSE

#ifdef SSE
#define ALLOC_EXPA true
#else
#define ALLOC_EXPA false
#endif

// Macros to avoid ifdef _OPENMP all over the code
#ifdef _OPENMP
#define CRITICAL _Pragma("omp critical")
#define PARFOR _Pragma("omp for") for
#define PARALLEL _Pragma("omp parallel num_threads(NUM_THREADS)")
#define NOPARALLEL _Pragma("omp parallel num_threads(1)")
#else
#define CRITICAL
#define PARFOR for
#define PARALLEL
#define NOPARALLEL
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int N = mxGetM(XVEC);
    double* restrict x = mxGetPr(XVEC);
    double* restrict n = mxGetPr(NVEC);
    double* restrict f = mxGetPr(FVEC);

    // pack parameters
    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    if(VERBOSE)
	mexPrintf("[Stresslet%s FG(G)] N=%d, P=%d\n", PER_STR,N,params.P);
    
    // allocate output array (M1,M2,M3,9)
    mwSize dims[4];
    dims[0] = params.dims[0];
    dims[1] = params.dims[1];
    dims[2] = params.dims[2];
    dims[3] = 9;
    H_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    size_t block = dims[0]*dims[1]*dims[2];
    double* H_base = mxGetPr(H_OUT);

    PARALLEL
    {
	// scratch arrays
	SE_FGG_work work;
	double* restrict q;
	CRITICAL
	{
	    SE_FGG_allocate_workspace(&work, &params,true,ALLOC_EXPA);    
	    q = mxMalloc(N*sizeof(double));
	}    
	// precompute
	NOPARALLEL // Disable work-sharing in FGG routines
	{
	    SE_FGG_base_gaussian(&work, &params);
#ifdef SSE
	    const SE_state st = {.x = x,  .q = NULL};
	    SE_FGG_expand_all(&work, &st, &params);	    
#endif
	}
	// loop over 9 output components
	PARFOR (int ii=0; ii<9; ii++)
	{
	    // ii is linear index
	    // ii = i1+3*i2
	    int i1 = ii % 3;
	    int i2 = (ii-i1)/3;
	    // Create product q = n_i1 * f_i2
	    for (int i=0; i<N; i++)
		q[i] = n[i+i1*N]*f[i+i2*N];
	    const SE_state st = {.x = x,  .q = q};
	    // Set output and clear work areas
	    double* H_per = &H_base[ii*block];
	    SE_fp_set_zero(work.H, SE_prod3(params.npdims));
	    SE_fp_set_zero(H_per, SE_prod3(params.dims));
	    // Now do the work
	    NOPARALLEL // Disable work-sharing in FGG routines
	    {
#ifdef SSE
#ifdef __AVX__
        SE_FGG_grid_split_AVX_dispatch(&work, &st, &params);
#else
		SE_FGG_grid_split_SSE_dispatch(&work, &st, &params);
#endif
#else
		SE_FGG_grid(&work, &st, &params);
#endif		
	    }
#ifdef THREE_PERIODIC
	    SE_FGG_wrap_fcn(H_per, &work, &params);
#endif    

#ifdef TWO_PERIODIC
	    SE2P_FGG_wrap_fcn(H_per, &work, &params);
#endif    
	}
	// All done, free work areas
	CRITICAL
	{
	    mxFree(q);
	    SE_FGG_free_workspace(&work);
	}
    }
}
