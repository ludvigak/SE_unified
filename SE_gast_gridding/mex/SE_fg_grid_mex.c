#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] 
#define Q   prhs[1] 
#define OPT prhs[2]

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
    SE_FGG_allocate_workspace(&work, &params,true,false);
    
    // allocate output array
    H_OUT = mxCreateNumericArray(3, params.dims, mxDOUBLE_CLASS, mxREAL);
    double* H_per = mxGetPr(H_OUT);
    SE_fp_set_zero(H_per, SE_prod3(params.dims));

    // coordinates and charges
    const SE_state st = {.x = x,  .q = q};

    if(VERBOSE)
	mexPrintf("[SE%s FG(G)] N=%d, P=%d\n",PER_STR,N,params.P);

    // now do the work
    SE_FGG_base_gaussian(&work, &params);
    SE_FGG_grid(&work, &st, &params);

#ifdef THREE_PERIODIC
    SE_FGG_wrap_fcn(H_per, &work, &params);
#endif    

#ifdef TWO_PERIODIC
    SE2P_FGG_wrap_fcn(H_per, &work, &params);
#endif    
    
    // done
    SE_FGG_free_workspace(&work);
}
