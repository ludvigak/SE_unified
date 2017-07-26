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
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>


#define pi     3.14159265358979323846264338327950288
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a > b) ? b : a)
#define EULER  0.577215664901532860606512090082402431042
#define EPSILON __DBL_EPSILON__
#define STACK_SIZE 100000  /*FIXME: This must be adpative.*/

double IncompBesselK0_Simpson(double tol, int *cnt, double bes_a, double bes_b, int der);
double computeK0(double a, double b);
double computeINCBK0(double a, double b, int der);

typedef struct {
  double a;
  double b;
} gsl_params;
double bessel_f(double x, void * p);
double bessel_f_der(double x, void * p);
double call_gsl_bessel_integrator(double a, double b, 
				  gsl_integration_workspace *w, int der);

void inline set_zero_4(double *x)
{ 
  x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0;
}

#endif
