#ifndef SE_DIRECT_H
#define SE_DIRECT_H

#define PI 3.141592653589793
#include <omp.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC mxMalloc
#define __PRINTF mexPrintf
#define __FREE mxFree
#else
#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#endif


typedef struct 
{
    double box[3];
    double xi;
    int layers; 
    double rc; 
} ewald_opts;

#endif
