#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define OPT prhs[0] 

#define Z plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    // pack parameters
    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, 0);

    // allocate output array
    int dims[3]={params.P, params.P, params.P};
    Z = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    SE_FGG_work work = {.H = NULL,  .zs = mxGetPr(Z)};

    if(VERBOSE)
	mexPrintf("[SE%s FG(S)] P=%d\n",PER_STR,params.P);

    // now do the work
    SE_FGG_base_gaussian(&work, &params);
    
    // done
}
