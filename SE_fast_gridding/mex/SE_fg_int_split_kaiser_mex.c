#include "mex.h"
#include "../SE_fgg.h"
#include "../SE_fkg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] // this arg is unused
#define HH  prhs[1] 
#define OPT prhs[2] 
#define ZX  prhs[3]
#define ZY  prhs[4]
#define ZZ  prhs[5]
#define IDX prhs[6]

#define PHI_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{

    const int N = mxGetM(IDX);
    const double* H_per = mxGetPr(HH);

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

    // output vector
    PHI_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    double* phi = mxGetPr(PHI_OUT);

    if(VERBOSE)
	mexPrintf("[SE%s FG(i)] N=%d, P=%d\n",PER_STR,N,params.P);

    // now do the work
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
#ifdef ONE_PERIODIC
      SE1P_FGG_extend_fcn(&work, H_per, &params);
#endif   
    SE_FKG_int_split_AVX_dispatch(phi, &work, &params);
    }

    
    // done
    SE_FGG_free_workspace(&work);
}
