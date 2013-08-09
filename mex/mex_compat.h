#ifndef __MEX_COMPAT_H__
#define __MEX_COMPAT_H__

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC mxMalloc
#define __REALLOC mxRealloc
#define __FREE mxFree
#define __PRINTF mexPrintf
#else
#include <stdlib.h>
#include <stdio.h>
#define __MALLOC malloc
#define __REALLOC realloc
#define __FREE free
#define __PRINTF printf
#endif

#endif
