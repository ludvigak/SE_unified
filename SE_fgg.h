#ifndef __SE_FGG_H
#define __SE_FGG_H

// System includes
#include "math.h"
#include <sys/time.h>
#include "emmintrin.h"

// Consttants and indexing 
#define PI 3.141592653589793
#define FGG_INF 1.79769e+308

#define __IDX3_CMAJ(II,IJ,IK,N1,N2) ( (II)+(IJ)*(N1)+(IK)*(N1)*(N2) )
#define __IDX3_RMAJ(II,IJ,IK,N2,N3) ( (II)*(N2)*(N3)+(IJ)*(N3)+(IK) )

// Select periodicty: must give -D<...> to compiler
#ifdef THREE_PERIODIC
#define __FGG_EXPA fgg_expansion_3p
#define PER_STR "3P"
#endif

#ifdef TWO_PERIODIC
#define __FGG_EXPA fgg_expansion_2p
#define PER_STR "2P"
#endif

// Maximal amount of Gaussian support (defined to help the compiler)
#define P_MAX 32

// Specific includes and defines for MEX-file compilation
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC malloc
#define __PRINTF mexPrintf
#define __FREE free
#define __IDX __IDX3_CMAJ
#else
#include <stdlib.h>
#include <stdio.h>
#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#define __IDX __IDX3_RMAJ
#endif

// display debug messages in the SSE dispatcher
#ifdef VERBOSE
#define __DISPATCHER_MSG(s) __PRINTF(s)
#else
#define __DISPATCHER_MSG(s) {}
#endif

inline int is_odd(int p)
// test if interger is odd or even
{
    return p&1;
}

// Return half of gaussian width:
//    q = (p-1)/2 if p is odd
//    q =  p/2    if p is even
inline int half(int p)
{
    return (is_odd(p) ? (p-1)/2 : p/2);
}

// Temporary, or auxillary arrays
typedef struct
{
    double* H;
    double* zs;

    double* zx;
    double* zy;
    double* zz;
    int* idx;

    int free_zs;
    int free_fgg_expa;

} SE_FGG_work;

// FGG parameters
typedef struct
{
    int N;
    int P;
    int P_half;
    int dims[3];
    int npdims[3];
    double c;
    double d;
    double h;
    double a;

} SE_FGG_params;

typedef struct
{
    int idx_on_grid, idx_in_array; 
} idx_reorder_t;

// Particle positions and charges
typedef struct
{
    double* x;
    double* q;

} SE_state;

// Fill parameter struct
void SE_FGG_pack_params(SE_FGG_params*, int, int, int, int, int, 
			double, double);
void SE2P_FGG_pack_params(SE_FGG_params*, int, int, int, int, int, 
			  double, double, double);

// Allocate workspace (malloc)
void SE_FGG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int, int);
double* SE_FGG_allocate_grid(const SE_FGG_params*);
double* SE_FGG_allocate_vec(int);

// Free workspace (free)
void SE_FGG_free_workspace(SE_FGG_work*);

// Particles to grid
void SE_FGG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_u8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P16(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Compute all FGG expansion vectors
void SE_FGG_expand_all(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Grid to particles
void SE_FGG_int(double*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_dispatch(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_u8(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P16(double*, const SE_FGG_work*, const SE_FGG_params*);

// Static Gaussian on P^3-grid
void SE_FGG_base_gaussian(SE_FGG_work*, const SE_FGG_params*);

// Wrap function to produce periodicity
void SE_FGG_wrap_fcn(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE2P_FGG_wrap_fcn(double*, const SE_FGG_work*, const SE_FGG_params*);

// Extend periodic function
void SE_FGG_extend_fcn(SE_FGG_work*, const double*, const SE_FGG_params*);
void SE2P_FGG_extend_fcn(SE_FGG_work*, const double*, const SE_FGG_params*);

// Randomize positions and charges (malloc)
void SE_init_system(SE_state*, const SE_FGG_params*);

// Free particles and charges (free)
void SE_free_system(SE_state*);

// Retrun time in seconds
double SE_gettime(void);

// Return product of elements in integer triplet
int SE_prod3(const int[3]);

// Set first N elements of array to floating-point zero
void SE_fp_set_zero(double*, int);

// Reorder particles according to closest grid point
void SE_FGG_reorder_system(SE_state*, const SE_FGG_work*, const SE_FGG_params*);

#endif
