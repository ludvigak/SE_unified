#ifndef __SE_FGG_H
#define __SE_FGG_H

// System includes
#include "math.h"
#include <sys/time.h>
#include "x86intrin.h"

// Consttants and indexing 
#define PI 3.14159265358979323846
#define FGG_INF 1.79769e+308

#define __IDX3_CMAJ(II,IJ,IK,N1,N2) ( (II)+(IJ)*(N1)+(IK)*(N1)*(N2) )
#define __IDX3_RMAJ(II,IJ,IK,N2,N3) ( (II)*(N2)*(N3)+(IJ)*(N3)+(IK) )

// Select periodicty: must give -D<...> to compiler
#ifdef THREE_PERIODIC
#define __FGG_EXPA fgg_expansion_3p
#define __FGG_EXPA_FORCE fgg_expansion_3p_force
#define __FGG_INDEX fgg_index_3p
#define PER_STR "3P"
#endif

#ifdef TWO_PERIODIC
#define __FGG_EXPA fgg_expansion_2p
#define __FGG_INDEX fgg_index_2p
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

// Malloc with 32 or 16-byte alignment (from intrinsics library)
#ifdef __AVX__
#define MEM_ALIGNED __attribute__((aligned(32)))
#define SE_FGG_MALLOC(sz) _mm_malloc((sz),32)
#define SE_FGG_FREE(sz) _mm_free((sz))
#else
#define MEM_ALIGNED __attribute__((aligned(16)))
#define SE_FGG_MALLOC(sz) _mm_malloc((sz),16)
#define SE_FGG_FREE(sz) _mm_free((sz))
#endif

// Print compile-time messages about which kernels will be used
#ifdef __AVX__
#if defined(__FMA__) && defined(__AVX2__)
#define AVX_FMA
//#pragma message ("Compiling vectorized kernels with AVX and FMA instructions.")
#else
//#pragma message ("Compiling vectorized kernels with AVX instructions (no FMA).")
#endif
#else
//#pragma message ("Compiling vectorized kernels with SSE instructions.")
#endif

// display debug messages in the SSE dispatcher
// only master thread prints when threaded
#ifdef VERBOSE
#ifdef _OPENMP
#define __DISPATCHER_MSG(s) _Pragma("omp master") \
    __PRINTF(s)
#else
#define __DISPATCHER_MSG(s) __PRINTF(s)
#endif
#else
#define __DISPATCHER_MSG(s) {}
#endif

static inline int is_odd(int p)
// test if interger is odd or even
{
    return p&1;
}
static inline int isnot_div_by_4(int p)
{
   return p&3;
}

// Return half of gaussian width:
//    q = (p-1)/2 if p is odd
//    q =  p/2    if p is even
static inline int half(int p)
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
  
  // extra space to compute the force
  double* zfx;
  double* zfy;
  double* zfz;
  
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
    double beta;

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
    double* phi;

} SE_state;

void SE_FGG_grid_kaiser(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Fill parameter struct
void SE_FGG_pack_params(SE_FGG_params*, int, int, int, int, int, 
			double, double);
void SE2P_FGG_pack_params(SE_FGG_params*, int, int, int, int, int, 
			  double, double, double);

// Allocate workspace (malloc)
void SE_FGG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int, int);
void SE_FGG_allocate_workspace_SSE_force(SE_FGG_work*, const SE_FGG_params*, int, int);
double* SE_FGG_allocate_grid(const SE_FGG_params*);
double* SE_FGG_allocate_vec(int);

// Free workspace (free)
void SE_FGG_free_workspace(SE_FGG_work*);
void SE_FGG_free_workspace_SSE_force(SE_FGG_work*);

// Particles to grid Potential
void SE_FGG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_u8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P16(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
#ifdef __AVX__
void SE_FGG_grid_split_AVX_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P16(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_u8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
#endif

// particles to grid force
void SE_FGG_grid_split_SSE_dispatch_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_u8_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P16_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
#ifdef __AVX__
void SE_FGG_grid_split_AVX_dispatch_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P16_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P8_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_u8_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
#endif

// Compute all FGG expansion vectors
void SE_FGG_expand_all(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_expand_all_SSE_force(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Grid to particles Potential
void SE_FGG_int(double*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_int_force(double*, const SE_FGG_work*, SE_state*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_dispatch(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_u8(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P16(double*, const SE_FGG_work*, const SE_FGG_params*);

// Grid to particles Force
#ifdef __AVX__
void SE_FGG_int_split_AVX_dispatch(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P8(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P16(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_u8(double*, const SE_FGG_work*, const SE_FGG_params*);
#endif

// Grid to particles Potential
void SE_FGG_int_split_SSE_dispatch_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_u8_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P16_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
#ifdef __AVX__
void SE_FGG_int_split_AVX_dispatch_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P16_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P8_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_u8_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
#endif

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
void SE_init_unit_system(SE_state*, const SE_FGG_params*);

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

// calculate energy
double calc_energy(SE_state, int);

#endif
