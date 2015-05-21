#include "math.h"
#include <sys/time.h>
#include "emmintrin.h"
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include "SE_fgg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef VERBOSE
#define __DISPATCHER_MSG(s) __PRINTF(s)
#else
#define __DISPATCHER_MSG(s) {}
#define VERBOSE 0
#endif

#define ONE_4PI_EPS0  138.9354491512000038

typedef double t_complex[2];

// input options
typedef struct 
{
  double m,c,box[3],h,xi,w,eta;
  int P,M,N;
} SE_opt;

// Unpacking params
void parse_params(SE_opt*, double);

// create k_space evctors
void k_vec(int, double*, double*, double*, double*);

// do the scaling
void scaling(double , double *, int , int, int, double*);


// products sr(scalar to real) rr (real to real) rc (real to complex)
// equivalent to .* in MATLAB. flag = 1 gives a.*b and flag = -1 gives a./b
void product_sr(double*, double, double*, int, int, int, int);
void product_rr(double*, double*, double*, int, int, int, int);
void product_rc(t_complex*, double*, t_complex*, int, int, int, int);


// Packing SE parameters
void
SE_FGG_FCN_params(SE_FGG_params*, const SE_opt*, int);

// calling gridding
void SE_fg_grid(double*, double*, int, SE_opt, double*);

// integration and interpolation
void SE_fgg_int(double*, double*, int, SE_opt, double*);

// integration and interpolation and calculate forces
void SE_fgg_int_force(double*, double *, double*, int, SE_opt, double*);

// 3d fft using fftw3 real to complex
void do_fft_r2c_3d(double*, t_complex*, int, int, int);
// 3d fft using fftw3 complex to real
void do_fft_c2r_3d(t_complex*, double*, int, int, int);
// 3d forward fft using fftw3 complex to complex
void do_fft_c2c_forward_3d(t_complex*, t_complex*, int, int, int);
// 3d backward fft using fftw3 complex to complex
void do_fft_c2c_backward_3d(t_complex*, t_complex*, int, int, int);



// printing results
void print_r1d(char*, double*, int, double,int);
void print_c1d(char*, t_complex*, int, double,int);
void print_r3d(char*, double*, int, int, int, double,int);
void print_c3d(char*, t_complex*, int, int, int, double,int);

void copy_r2c(double*, t_complex*, int);
void copy_c2r(t_complex*, double*, int);

double
SE_init_params(int*, int*, double*, double, double,
               double, const double*, int);
double LambertW(const double);
