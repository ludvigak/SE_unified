#include "math.h"
#include "string.h"
#include "time.h"

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

#define PI 3.141592653589793

#ifndef VERBOSE
#define VERBOSE 0
#endif

void  get_rs_triplets (const double* restrict x, const double* restrict nvec, int N,
		       const double* restrict box, double xi, double rc, int nlhs,
		       int* restrict *row_p, int* restrict *col_p, double* restrict val[3][3],
		       int* restrict *buck_size_p, int* restrict *idx_in_array_p, int* numel_p
		       );
