#ifndef __STRESSLET_REAL_RC_H__
#define __STRESSLET_REAL_RC_H__

#include "math.h"
#include "string.h"
#include "time.h"
#include "sys/time.h"
#include "unistd.h"

#include "mex_compat.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif


#define DELTA(tic,toc) ((toc.tv_sec  - tic.tv_sec) * 1000000u + toc.tv_usec - tic.tv_usec) / 1.e6

#define PI 3.141592653589793

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

#endif

