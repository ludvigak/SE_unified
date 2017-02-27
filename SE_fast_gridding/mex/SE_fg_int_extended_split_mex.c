#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] // this arg is unused
#define HH  prhs[1] 
#define OPT prhs[2] 
#define ZS  prhs[3]
#define ZX  prhs[4]
#define ZY  prhs[5]
#define ZZ  prhs[6]
#define IDX prhs[7]

#define PHI_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{

    const int N = mxGetM(IDX);

    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    SE_FGG_work work;
    // attach pre-computed quantities
    work.zs = mxGetPr(ZS);
    work.zx = mxGetPr(ZX);
    work.zy = mxGetPr(ZY);
    work.zz = mxGetPr(ZZ);
    work.idx = (int*)mxGetData(IDX);
    work.H = mxGetPr(HH);

    // output vector
    PHI_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    double* phi = mxGetPr(PHI_OUT);

    if(VERBOSE)
	mexPrintf("[SE%s FG(i)] N=%d, P=%d\n",PER_STR,N,params.P);

    if(N==1)
    {
	// Don't thread for single target
#ifdef __AVX__	
	SE_FGG_int_split_AVX_dispatch(phi, &work, &params);
#else
	SE_FGG_int_split_SSE_dispatch(phi, &work, &params);
#endif
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
	{
	    // now do the work
#ifdef __AVX__	
	    SE_FGG_int_split_AVX_dispatch(phi, &work, &params);
#else
	    SE_FGG_int_split_SSE_dispatch(phi, &work, &params);
#endif
	
	}
    }
}
