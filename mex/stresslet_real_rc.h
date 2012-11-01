#include "math.h"
#include "string.h"
#include "time.h"
#include "sys/time.h"
#include "unistd.h"

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

#ifdef _OPENMP
#include "omp.h"
#endif


#define DELTA(tic,toc) ((toc.tv_sec  - tic.tv_sec) * 1000000u + toc.tv_usec - tic.tv_usec) / 1.e6

#define PI 3.141592653589793

#ifndef VERBOSE
#define VERBOSE 0
#endif

void  get_rs_triplets (const double* restrict x, const double* restrict nvec, int N,
		       const double* restrict box, double xi, double rc, int nlhs,
		       int* restrict *row_p, int* restrict *col_p, double* restrict val[3][3],
		       int* restrict *buck_size_p, int* restrict *idx_in_array_p, int* numel_p
		       );

void  compute_rsrc_direct (const double* restrict x, 
			   const double* restrict nvec, 
			   const double* restrict qvec, 
			   int N,
			   const double* restrict box, 
			   double xi, 
			   double rc, 
			   double* restrict *phi_p
			   );

