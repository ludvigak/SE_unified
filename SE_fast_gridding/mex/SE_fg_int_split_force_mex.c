#include "mex.h"
#include "../SE_fgg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0]
#define HH  prhs[1] 
#define OPT prhs[2] 
#define ZS  prhs[3]
#define ZX  prhs[4]
#define ZY  prhs[5]
#define ZZ  prhs[6]
#define ZFX prhs[7]
#define ZFY prhs[8]
#define ZFZ prhs[9]
#define IDX prhs[10]

#define FORCE_OUT plhs[0]  // Output

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
  SE_FGG_allocate_workspace(&work, &params, false, false);

  // attach pre-computed quantities
  work.zs = mxGetPr(ZS);
  work.zx = mxGetPr(ZX);
  work.zy = mxGetPr(ZY);
  work.zz = mxGetPr(ZZ);
  work.zfx = mxGetPr(ZFX);
  work.zfy = mxGetPr(ZFY);
  work.zfz = mxGetPr(ZFZ);
  work.idx = (int*)mxGetData(IDX);

  // output vector
  FORCE_OUT = mxCreateDoubleMatrix(N,3,mxREAL);
  double* force = mxGetPr(FORCE_OUT);

  // coordinates and charges
  SE_state st = {.x = x, .q = NULL};
  
  if(VERBOSE)
    mexPrintf("[SE%s FG(i)] N=%d, P=%d\n",PER_STR,N,params.P);

  // now do the work
#ifdef THREE_PERIODIC
  SE_FGG_extend_fcn(&work, H_per, &params);
#endif

#ifdef TWO_PERIODIC
  SE2P_FGG_extend_fcn(&work, H_per, &params);
#endif

#ifdef ONE_PERIODIC
  SE1P_FGG_extend_fcn(&work, H_per, &params);
#endif
 
#ifdef __AVX__
  SE_FGG_int_split_AVX_dispatch_force(force, &st, &work, &params);
#elif __SSE4_2__
  SE_FGG_int_split_SSE_dispatch_force(force, &st, &work, &params);
#endif

  // done
  SE_FGG_free_workspace(&work);
}
