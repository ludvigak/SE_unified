#ifndef __MEX_COMPAT_H__
#define __MEX_COMPAT_H__

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC mxMalloc
#define __CALLOC mxCalloc
#define __REALLOC mxRealloc
#define __FREE mxFree
#define __PRINTF mexPrintf
#define __ERROR(msg) { \
	char ebuffer[1024];						\
	snprintf(ebuffer, 1024, "%s at %s, line %d.\n", msg, __FILE__, __LINE__); \
	mexErrMsgTxt(ebuffer);						\
    }
#define __WARNING(msg) { \
	char wbuffer[1024];						\
	snprintf(wbuffer, 1024, "%s at %s, line %d.\n", msg, __FILE__, __LINE__); \
	mexWarnMsgTxt(wbuffer);						\
    }
#define ASSERT(expr, err)						\
    if(! (expr))							\
	__ERROR("Assertion failed: (" # expr ") " # err);	
#else
#include <stdlib.h>
#include <stdio.h>
#define __MALLOC malloc
#define __CALLOC calloc
#define __REALLOC realloc
#define __FREE free
#define __PRINTF printf
#define __ERROR(msg) assert(0 && msg)
#define __WARNING(msg) printf("[WARNING] %s at %s, line %d.\n", msg, __FILE__, __LINE__);
#define ASSERT(expr, err) 						\
	assert( (expre) && "Assert (" # expr ") failed: " # err "\n");	
#endif

#endif
