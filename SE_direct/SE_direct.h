#ifndef SE_DIRECT_H
#define SE_DIRECT_H

#define PI 3.141592653589793

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

void unpack_opt(ewald_opts* opt, const mxArray* mx_opt);
void index_translation(int* idx, const double* idx_d, int N);
void SE2P_direct_real(double*, const int*, int,
		      const double*, const double*, int, 
		      const ewald_opts);

void SE2P_direct_real_rc(double*, const int*, int,
			 const double*, const double*, int, 
			 const ewald_opts);

void SE2P_direct_fd(double*, const int*, int,
		    const double *, const double*, int, 
		    const ewald_opts);

void SE2P_direct_k0(double*, const int*, int,
		    const double *, const double*, int, 
		    const ewald_opts);

void SE2P_direct_self(double*, const int*, int,
		      const double*, int, const ewald_opts);

void SE3P_direct_real(double*, const int*, int,
		      const double*, const double*, int, 
		      const ewald_opts);

void SE3P_direct_real_rc(double*, const int*, int,
			 const double*, const double*, int, 
			 const ewald_opts);

void SE3P_direct_fd(double*, const int*, int,
		    const double *, const double*, int, 
		    const ewald_opts);

void SE3P_direct_self(double*, const int*, int,
		      const double*, int, const ewald_opts);

#endif
