#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define HIN prhs[0] 
#define OPT prhs[1] 

#define HOUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int N = 1;
    const double* H_per = mxGetPr(HIN);

    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);
    size_t dims[3] = {params.npdims[0], params.npdims[1], params.npdims[2]};
    HOUT = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    SE_FGG_work work;
    work.H = mxGetPr(HOUT);

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
#ifdef THREE_PERIODIC
	SE_FGG_extend_fcn(&work, H_per, &params);
#endif
	
#ifdef TWO_PERIODIC
	SE2P_FGG_extend_fcn(&work, H_per, &params);
#endif
    }
}
