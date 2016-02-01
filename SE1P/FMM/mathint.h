#ifndef __MATHINT_H_
#define __MATHINT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include <stdbool.h>
#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif


#define pi     3.14159265358979323846264338327950288
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a > b) ? b : a)
#define EPS    1.110223024625157e-16
#define EULER  0.577215664901532860606512090082402431042
#define GSL_LOG_DBL_MIN 7.083964185322641e+02


void glwt_fast(int, double, double, double*, double*);
void glwt(int, double, double, double*, double*);
void trapz_wt(int, double*, double*);
double expint(int, double);
double matlab_expint(double);
double expint_opt(double);
double expint_log_euler(double);
double IncompBesselK0_int_inf(double a, double b, int n);
double IncompBesselK0_int(double a, double b, int n, double* x, double* w);


struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
};
typedef struct cheb_series_struct cheb_series;

double gsl_expint_log_euler_1(  const double*, const double*, const int);
double gsl_expint_log_euler_4(  const double*, const double*, const int);
double gsl_expint_log_euler_32( const double*, const double*, const int);
double gsl_expint_log_euler_MAX(const double*, const double*, const int);
void constants(double *, double *, double *, double *, 
	       int, const double, const double, const int);
double k0_term(double *,double *,double *,double *,
	       const double * restrict , const double * restrict, int, int, 
	       const int, const double, const double);
void inline set_zero_4(double *x)
{ 
  x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0;
}

#endif
