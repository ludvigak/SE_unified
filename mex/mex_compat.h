#ifndef __MEX_COMPAT_H__
#define __MEX_COMPAT_H__

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC mxMalloc
#define __CALLOC mxCalloc
#define __REALLOC mxRealloc
#define __FREE mxFree
#define __PRINTF mexPrintf
#define __ERROR mexErrMsgTxt
#define ASSERT(expr, err)						\
    if(! (expr))							\
	mexErrMsgTxt("Assertion failed: (" # expr ") " # err "\n");	
#else
#include <stdlib.h>
#include <stdio.h>
#define __MALLOC malloc
#define __CALLOC calloc
#define __REALLOC realloc
#define __FREE free
#define __PRINTF printf
#define __ERROR(msg) assert(0 && msg)
#define ASSERT(expr, err) 						\
	assert( (expre) && "Assert (" # expr ") failed: " # err "\n");	
#endif

#endif
