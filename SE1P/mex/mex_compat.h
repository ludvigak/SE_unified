#ifndef __MEX_COMPAT_H__
#define __MEX_COMPAT_H__

#define GiB *((size_t) 1024*1024*1024)

/** Maximum total amount of memory to malloc in applications,
    not strictly followed, but tries to avoid some swap-deaths. */
#define MALLOC_MAX 4 GiB


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
#define __FLUSH() mexEvalString("drawnow;")
#else
#include <stdlib.h>
#include <stdio.h>
#include "assert.h"
#define __MALLOC malloc
#define __CALLOC calloc
#define __REALLOC realloc
#define __FREE free
#define __PRINTF printf
#define __ERROR(msg) assert(0 && msg)
#define __WARNING(msg) printf("[WARNING] %s at %s, line %d.\n", msg, __FILE__, __LINE__);
#define ASSERT(expr, err) 						\
	assert( (expr) && "Assert (" # expr ") failed: " # err "\n");	
#define __FLUSH() 
#endif

#endif
