#include "mex.h"
#include "../SE_fgg.h"
#include "../SE_fkg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] 
#define Q   prhs[1] 
#define OPT prhs[2]
#define ZX  prhs[3]
#define ZY  prhs[4]
#define ZZ  prhs[5]
#define IDX prhs[6]

#define H_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int N = mxGetM(X);
    double* restrict x = mxGetPr(X);
    double* restrict q = mxGetPr(Q);

    // pack parameters
    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    // scratch arrays
    SE_FGG_work work;
    SE_FKG_allocate_workspace(&work, &params, false);

    // attach pre-computed quantities
    work.zx = mxGetPr(ZX);
    work.zy = mxGetPr(ZY);
    work.zz = mxGetPr(ZZ);
    work.idx = (int*)mxGetData(IDX);
    
    // allocate output array
    size_t dims[3] = {params.dims[0], params.dims[1], params.dims[2]};
    H_OUT = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* H_per = mxGetPr(H_OUT);
    SE_fp_set_zero(H_per, SE_prod3(params.dims));

    // coordinates and charges
    const SE_state st = {.x = x,  .q = q};

    if(VERBOSE)
	mexPrintf("[SE%s FG(g)] N=%d, P=%d\n",PER_STR,N,params.P);

    // now do the work
    /* SE_FKG_grid_split_SSE_dispatch(&work, &st, &params); */
    SE_FKG_grid_split_AVX_dispatch(&work, &st, &params);
        
#ifdef THREE_PERIODIC
    SE_FGG_wrap_fcn(H_per, &work, &params);
#endif    



    // done
    SE_FGG_free_workspace(&work);
}
