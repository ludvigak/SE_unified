#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] 
#define HH  prhs[1] 
#define OPT prhs[2] 

#define PHI_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int N = mxGetM(X);
    double* restrict x = mxGetPr(X);
    const double* H_per = mxGetPr(HH);

    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    // scratch arrays
    SE_FGG_work work;
    SE_FGG_allocate_workspace(&work, &params,true,false);

    // output vector
    PHI_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    double* phi = mxGetPr(PHI_OUT);

    // coordinates and charges
    const SE_state st = {.x = x,  .q = NULL};

    if(VERBOSE)
	mexPrintf("[SE%s FG(I)] N=%d, P=%d\n",PER_STR,N,params.P);

    // now do the work
    SE_FGG_base_gaussian(&work, &params);

#ifdef THREE_PERIODIC
    SE_FGG_extend_fcn(&work, H_per, &params);
#endif

#ifdef TWO_PERIODIC
    SE2P_FGG_extend_fcn(&work, H_per, &params);
#endif

    SE_FGG_int(phi, &work, &st, &params);

    // done
    SE_FGG_free_workspace(&work);
}
