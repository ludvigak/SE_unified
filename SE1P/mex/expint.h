#ifndef EXPINT_H
#define EXPINT_H
/*------------------------------------------------------------------------
 *This file contains different implementations of the expint and ein functions
 *based on the MATLAB implementation of the expint function and gsl gnu library.
 *The SSE and AVX implementations of both function are available. For double and 
 *single precisions, compile with -DDOUBLE and -DSINGLE.
 *
 *Usage:
 * expint(x), ein(x), gsl_expint(x), gsl_ein(x)
 * _mm_expint_ps(X),     _mm_expint_pd(X),    _mm_ein_ps(X),    _mm_ein_pd(X)
 * _mm256_expint_ps(X),  _mm256_expint_pd(X), _mm256_ein_ps(X), _mm256_ein_pd(X)
 *
 * Davoud Saffar, davoudss@kth.se, Apr. 2016.
 */

#include "stdbool.h"
#include "float.h"
#include "math_x86.h"

/*Define alignment*/
#ifdef __AVX__
#define MEM_ALIGN __attribute__((aligned(32)));
#else
#define MEM_ALIGN __attribute__((aligned(16)));
#endif

#ifdef DOUBLE
typedef double real;
#define EPS     __DBL_EPSILON__
#else
typedef float real;
#define EPS     __FLT_EPSILON__
#endif

#if defined __GNUC__ && __AVX__
#ifdef DOUBLE
#define _mm256_log_pd se_mm256_log_pd
#define _mm256_exp_pd se_mm256_exp_pd
#else
#define _mm256_log_ps se_mm256_log_ps
#define _mm256_exp_ps se_mm256_exp_ps
#endif
#endif

#define EULER   0.577215664901532860606512090082402431042
#define INFTY   INFINITY



static double E1_data[16] = {
  -0.03739021479220279500,
   0.04272398606220957700,
  -0.13031820798497005440,
   0.01441912402469889073,
  -0.00134617078051068022,
   0.00010731029253063780,
  -0.00000742999951611943,
   0.00000045377325690753,
  -0.00000002476417211390,
   0.00000000122076581374,
  -0.00000000005485141480,
   0.00000000000226362142,
  -0.00000000000008635897,
   0.00000000000000306291,
  -0.00000000000000010148,
   0.00000000000000000315
};

static double E4_data[25] = {
  -0.605773246640603460,
  -0.112535243483660900,
   0.013432266247902779,
  -0.001926845187381145,
   0.000309118337720603,
  -0.000053564132129618,
   0.000009827812880247,
  -0.000001885368984916,
   0.000000374943193568,
  -0.000000076823455870,
   0.000000016143270567,
  -0.000000003466802211,
   0.000000000758754209,
  -0.000000000168864333,
   0.000000000038145706,
  -0.000000000008733026,
   0.000000000002023672,
  -0.000000000000474132,
   0.000000000000112211,
  -0.000000000000026804,
   0.000000000000006457,
  -0.000000000000001568,
   0.000000000000000383,
  -0.000000000000000094,
   0.000000000000000023
};

static double E32_data[26] = {
  -0.18929180007530170,
  -0.08648117855259871,
   0.00722410154374659,
  -0.00080975594575573,
   0.00010999134432661,
  -0.00001717332998937,
   0.00000298562751447,
  -0.00000056596491457,
   0.00000011526808397,
  -0.00000002495030440,
   0.00000000569232420,
  -0.00000000135995766,
   0.00000000033846628,
  -0.00000000008737853,
   0.00000000002331588,
  -0.00000000000641148,
   0.00000000000181224,
  -0.00000000000052538,
   0.00000000000015592,
  -0.00000000000004729,
   0.00000000000001463,
  -0.00000000000000461,
   0.00000000000000148,
  -0.00000000000000048,
   0.00000000000000016,
  -0.00000000000000005
};

// valid for 0<x<=1
real
cheb_eval_e1(const real x)
{
  int j;
  real d  = 0.0;
  real dd = 0.0;
  // this is cancelled as a result of adding log term
  real ln_term = -log(x);
  real y2 = 2.0 * x;

#ifdef DOUBLE
  for(j = 15; j>=1; j--) {
#else
  for(j = 11; j>=1; j--) {
#endif
    real temp = d;
    d = y2*d - dd + E1_data[j];
    dd = temp;
  }
      
  d = x*d - dd + 0.5 * E1_data[0];
      
  // add log term (cancelled with ln_term) and EULER constant
  return -0.6875+x+d+ln_term;
}

// valid for 1<x<=4
real
cheb_eval_e4(const real x)
{
  int j;
  real s;

  real d  = 0.0;
  real dd = 0.0;

  s = 1.0/x*exp(-x);
  real x0 = (8.0/x-5.0)/3.0;

  real y2 = 2.0 * x0;
#ifdef DOUBLE
  for(j = 24; j>=1; j--) 
#else
  for(j = 11; j>=1; j--) 
#endif
    {
      real temp = d;
      d = y2*d - dd + E4_data[j];
      dd = temp;
    }
      
  d = x0*d - dd + 0.5 * E4_data[0];

  // add log term and EULER constant
  return s*(1.0+d);
}

// valid for 4<x<=32
real
cheb_eval_e32(const real x)
{
  int j;
  real s;
  real d  = 0.0;
  real dd = 0.0;
  s = 1.0/x*exp(-x);
  real x0 = 8.0/x-1.0;
	
  real y2 = 2.0 * x0;
#ifdef DOUBLE
  for(j = 25; j>=1; j--) 
#else
  for(j = 11; j>=1; j--) 
#endif
    {
      real temp = d;
      d = y2*d - dd + E32_data[j];
      dd = temp;
    }
	
  d = x0*d - dd + 0.5 * E32_data[0];
  // add log term and EULER constant
  return s*(1.0+d);

}

real
cheb_eval_eMAX(real x)
{
  return 0.0;
}

real gsl_expint(const real x)
{
  if(x<=__DBL_EPSILON__)
    return 0;
  else if(x<=1)
    return cheb_eval_e1(x);
  else if(x<=4)
    return cheb_eval_e4(x);
  else if(x<=32)
    return cheb_eval_e32(x);
  else
    return cheb_eval_eMAX(x);
  return 0;
}


/* sum up 4 elements in a vector, used for SIMD op.*/
inline double sum4(double*p){
  /*Unexpected behaviour if p has more or less
   *elements than 4 */
  return (p[0]+p[1]+p[2]+p[3]);
}

inline float sum8(float*p){
  /*Unexpected behaviour if p has more or less
   *elements than 4 */
  return (p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7]);
}

double expint(double x)
{
  int j;
  double am1,am2,bm1,bm2,f,oldf,pterm,term,y,egamma,a,b,alpha,beta;
  double INF = INFINITY;

  if(x<3.0162691827353)
    {
      egamma=EULER;
      y = -egamma-log(x);
      j = 1;
      pterm = x;
      term = x;
      
      while(fabs(term)>EPS)
	{
	  y += term;
	  j ++;
	  pterm = -x*pterm/(double)j;
	  term  = pterm/(double) j;
	}
    }
  else
    {
      am2 = 0.;bm2 = 1.; 
      am1 = 1.;bm1 = x;
      
      f    = am1 / bm1;
      oldf = INF;
      j    = 2;

      while( fabs(f-oldf) > 100.*EPS*fabs(f))
	{
	  alpha = (double) j/2.;
	  
	  /* calculate A(j), B(j), and f(j)*/
	  a = am1 + alpha * am2;
	  b = bm1 + alpha * bm2;
   
	  /* save new normalized variables for next pass through the loop
	   *  note: normalization to avoid overflow or underflow*/
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1.0;
	  
	  f = am1;
	  j++;
   
	  /* calculate the coefficients for j odd*/
	  alpha = (double) (j-1)/2.;
	  beta = x;
	  a = beta  * am1 + alpha * am2;
	  b = beta  * bm1 + alpha * bm2;
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1.0;
	  oldf= f;
	  f   = am1;
	  j++;
	} /* end of while loop*/
      
      y = exp(-x)*f;

    }

  return y;
}

real ein(real x)
{
  real egamma = 0.57721566490153286061;
  if(fabs(x)<EPS)
    return 0.0;
  else if(fabs(x)>32)
    return log(x)+egamma;
  else
    return expint(x)+log(x)+egamma;
}


real gsl_ein(const real x)
{
  if(x<=__DBL_EPSILON__)
    return 0;
  else if(x<=1)
    return cheb_eval_e1(x)+log(x)+EULER;
  else if(x<=4)
    return cheb_eval_e4(x)+log(x)+EULER;
  else if(x<=32)
    return cheb_eval_e32(x)+log(x)+EULER;
  else
    return cheb_eval_eMAX(x)+log(x)+EULER;
  return 0;
}

#ifdef  __AVX__
#include "math_x86.h"
/*An extra test should be added for a more general 
 *routine to check whether any of the double precision
 *elements are zero or not. This can be done by comparing 
 *against EPS and setting zero values to EPS. e.g.,
 *Y = _mm256_cmp_pd(Z,EPS,_CMP_LE_OQ);
 *Z = _mm256_add_pd(_mm256_and_pd(ONE16,Y),Z);
 *We are safe since we do not use this routine for 
 *computing the direct interactions of the particles 
 *in the same box.*/
/*This routine can compute each pack of 4 double precision
 *values in 1e-7sec.*/
__m256d _mm256_expint_pd(__m256d Z)
{
  int MAXITER = 25;
  double branch_pred     = 3.0162691827353;
  __m256d GAMMA          = _mm256_set1_pd(EULER);
  __m256d C0_BRANCH_PRED = _mm256_set1_pd(branch_pred);
  __m256d ONE            = _mm256_set1_pd(1.0);
  __m256d ZERO           = _mm256_setzero_pd();
  __m256d ERR,Z0,Z1;
  
  /* branch prediction */
  __m256d C1 = _mm256_cmp_pd(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  double s[4] MEM_ALIGN;
  double f[4] MEM_ALIGN;
  
  
  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = true;
  int j = 1;
  /* If one of the double precision values is large, this loop
   * reaches the MAXITER. Therefore, we set all the values bigger than
   * C0_BRANCH_PRED, to this value. The same is done for the continued 
   * fraction iteration.
   */
  Z0 = _mm256_min_pd(Z,C0_BRANCH_PRED);
  
  __m256d Y     = _mm256_sub_pd(ZERO,_mm256_add_pd(GAMMA,_mm256_log_pd(Z0)));
  __m256d PTERM = Z0;
  __m256d TERM  = Z0;
  
  /* The output will be in Y.*/
  while (term==true && j<MAXITER){
    Y     = _mm256_add_pd(Y, TERM);
    j++;
    PTERM = _mm256_mul_pd(Z0,_mm256_mul_pd(PTERM,_mm256_set1_pd(-1.0/j)));
    TERM  = _mm256_mul_pd(PTERM,_mm256_set1_pd(1.0/j));
    
    _mm256_store_pd(s,TERM);
    term = (fabs(s[0])>EPS) || (fabs(s[1])>EPS) || (fabs(s[2])>EPS) || (fabs(s[3])>EPS);
    
  }
   
  _mm256_store_pd(s,_mm256_sub_pd(Z,Z0));
  term = (fabs(s[0])<EPS) && (fabs(s[1])<EPS) && (fabs(s[2])<EPS) && (fabs(s[3])<EPS);
  if(term==true)
    return Y;


  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  double alpha;
  bool err = 1;
  Z1 = _mm256_max_pd(Z,C0_BRANCH_PRED);

  __m256d MINUS= _mm256_set1_pd(-1.0);
  __m256d AM2 = ZERO;
  __m256d BM2 = ONE;
  __m256d AM1 = ONE;
  __m256d BM1 = Z1;
  __m256d F   = _mm256_div_pd(AM1,BM1);
  __m256d OLDF= _mm256_set1_pd(INFTY);
  __m256d ALPHA, A, B, INVB, BETA;

  while (err==1 && j<MAXITER+10){
    /* calculate the coefficients of the recursion formulas for j even*/
    alpha = j/2.0;
    ALPHA = _mm256_set1_pd(alpha);
    A = _mm256_add_pd(AM1, _mm256_mul_pd(ALPHA,AM2));
    B = _mm256_add_pd(BM1, _mm256_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_pd(B);
    AM2  = _mm256_mul_pd(AM1,INVB);
    BM2  = _mm256_mul_pd(BM1,INVB);
    AM1  = _mm256_mul_pd(A,INVB);
    BM1  = ONE;

    F    = AM1;
    j++;

    /*calculate the coefficients for j odd*/
    alpha = (j-1.0)/2.0;
    ALPHA = _mm256_set1_pd(alpha);
    BETA  = Z1;
    A = _mm256_add_pd(_mm256_mul_pd(BETA,AM1), _mm256_mul_pd(ALPHA,AM2));
    B = _mm256_add_pd(_mm256_mul_pd(BETA,BM1), _mm256_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_pd(B);
    AM2  = _mm256_mul_pd(AM1,INVB);
    BM2  = _mm256_mul_pd(BM1,INVB);
    AM1  = _mm256_mul_pd(A,INVB);
    BM1  = ONE;

    OLDF = F;
    F    = AM1;
    j++;

    ERR = _mm256_sub_pd(F,OLDF);
    _mm256_store_pd(s,ERR);
    _mm256_store_pd(f,F);
    err = (fabs(s[0]/f[0])>EPS) || (fabs(s[1]/f[1])>EPS) || 
	  (fabs(s[2]/f[2])>EPS) || (fabs(s[3]/f[3])>EPS);

  }

  /* The output is in F*/
  F = _mm256_mul_pd(F, _mm256_exp_pd(_mm256_mul_pd(MINUS,Z1)));

  /* Return the corresponding answer */
  Y = _mm256_or_pd( _mm256_and_pd(C1,Y),_mm256_andnot_pd(C1,F) );

  return Y;
}


__m256d _mm256_ein_pd(__m256d Z)
{
  int MAXITER = 25;
  double branch_pred     = 3.0162691827353;
  __m256d GAMMA          = _mm256_set1_pd(EULER);
  __m256d C0_BRANCH_PRED = _mm256_set1_pd(branch_pred);
  __m256d ONE            = _mm256_set1_pd(1.0);
  __m256d ZERO           = _mm256_setzero_pd();
  __m256d LOG_P_EULER    = _mm256_add_pd(_mm256_log_pd(Z),GAMMA);
  __m256d ERR,Z0,Z1;
  
  /* branch prediction */
  __m256d C1 = _mm256_cmp_pd(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  double s[4] MEM_ALIGN;
  double f[4] MEM_ALIGN;
  
  
  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = true;
  int j = 1;
  /* If one of the double precision values is large, this loop
   * reaches the MAXITER. Therefore, we set all the values bigger than
   * C0_BRANCH_PRED, to this value. The same is done for the continued 
   * fraction iteration.
   */
  Z0 = _mm256_min_pd(Z,C0_BRANCH_PRED);
  
  __m256d Y = _mm256_sub_pd(ZERO,_mm256_add_pd(GAMMA,_mm256_log_pd(Z0)));
  __m256d PTERM = Z0;
  __m256d TERM  = Z0;
  
  /* The output will be in Y.*/
  while (term==true && j<MAXITER){
    Y     = _mm256_add_pd(Y, TERM);
    j++;
    PTERM = _mm256_mul_pd(Z0,_mm256_mul_pd(PTERM,_mm256_set1_pd(-1.0/j)));
    TERM  = _mm256_mul_pd(PTERM,_mm256_set1_pd(1.0/j));
    
    _mm256_store_pd(s,TERM);
    term = (fabs(s[0])>EPS) || (fabs(s[1])>EPS) || (fabs(s[2])>EPS) || (fabs(s[3])>EPS);
    
  }
  
  /* add the contribution of the log end euler constant */
  Y = _mm256_add_pd( LOG_P_EULER, Y);
  
  _mm256_store_pd(s,_mm256_sub_pd(Z,Z0));
  term = (fabs(s[0])<EPS) && (fabs(s[1])<EPS) && (fabs(s[2])<EPS) && (fabs(s[3])<EPS);
  if(term==true)
    return Y;


  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  double alpha;
  bool err = 1;
  Z1 = _mm256_max_pd(Z,C0_BRANCH_PRED);

  __m256d MINUS= _mm256_set1_pd(-1.0);
  __m256d AM2 = ZERO;
  __m256d BM2 = ONE;
  __m256d AM1 = ONE;
  __m256d BM1 = Z1;
  __m256d F   = _mm256_div_pd(AM1,BM1);
  __m256d OLDF= _mm256_set1_pd(INFTY);
  __m256d ALPHA, A, B, INVB, BETA;

  while (err==1 && j<MAXITER+10){
    /* calculate the coefficients of the recursion formulas for j even*/
    alpha = j/2.0;
    ALPHA = _mm256_set1_pd(alpha);
    A = _mm256_add_pd(AM1, _mm256_mul_pd(ALPHA,AM2));
    B = _mm256_add_pd(BM1, _mm256_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_pd(B);
    AM2  = _mm256_mul_pd(AM1,INVB);
    BM2  = _mm256_mul_pd(BM1,INVB);
    AM1  = _mm256_mul_pd(A,INVB);
    BM1  = ONE;

    F    = AM1;
    j++;

    /*calculate the coefficients for j odd*/
    alpha = (j-1.0)/2.0;
    ALPHA = _mm256_set1_pd(alpha);
    BETA  = Z1;
    A = _mm256_add_pd(_mm256_mul_pd(BETA,AM1), _mm256_mul_pd(ALPHA,AM2));
    B = _mm256_add_pd(_mm256_mul_pd(BETA,BM1), _mm256_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_pd(B);
    AM2  = _mm256_mul_pd(AM1,INVB);
    BM2  = _mm256_mul_pd(BM1,INVB);
    AM1  = _mm256_mul_pd(A,INVB);
    BM1  = ONE;

    OLDF = F;
    F    = AM1;
    j++;

    ERR = _mm256_sub_pd(F,OLDF);
    _mm256_store_pd(s,ERR);
    _mm256_store_pd(f,F);
    err = (fabs(s[0]/f[0])>EPS) || (fabs(s[1]/f[1])>EPS) || 
	  (fabs(s[2]/f[2])>EPS) || (fabs(s[3]/f[3])>EPS);

  }

  F = _mm256_mul_pd(F, _mm256_exp_pd(_mm256_mul_pd(MINUS,Z1)));

  /* The output is in F*/
  /* add the contribution of the log end euler constant */
  F = _mm256_add_pd( LOG_P_EULER, F);

  /* Return the corresponding answer */
  Y = _mm256_or_pd( _mm256_and_pd(C1,Y),_mm256_andnot_pd(C1,F) );

  return Y;
}

#ifdef SINGLE
// single precision routines
__m256 _mm256_expint_ps(__m256 Z)
{
  int MAXITER = 10;
  float branch_pred       = 3.0162691827353;
  __m256 GAMMA          = _mm256_set1_ps(EULER);
  __m256 C0_BRANCH_PRED = _mm256_set1_ps(branch_pred);
  __m256 ONE            = _mm256_set1_ps(1.0);
  __m256 ZERO           = _mm256_setzero_ps();
  __m256 ERR,Z0,Z1;
  
  /* branch prediction */
  __m256 C1 = _mm256_cmp_ps(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  float s[8] MEM_ALIGN;
  float f[8] MEM_ALIGN;
  
  
  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = true;
  int j = 1;
  /* If one of the single precision values is large, this loop
   * reaches the MAXITER. Therefore, we set all the values bigger than
   * C0_BRANCH_PRED, to this value. The same is done for the continued 
   * fraction iteration.
   */
  Z0 = _mm256_min_ps(Z,C0_BRANCH_PRED);
  
  __m256 Y     = _mm256_sub_ps(ZERO,_mm256_add_ps(GAMMA,_mm256_log_ps(Z0)));
  __m256 PTERM = Z0;
  __m256 TERM  = Z0;
  
  /* The output will be in Y.*/
  while (term==true && j<MAXITER){
    Y     = _mm256_add_ps(Y, TERM);
    j++;
    PTERM = _mm256_mul_ps(Z0,_mm256_mul_ps(PTERM,_mm256_set1_ps(-1.0/j)));
    TERM  = _mm256_mul_ps(PTERM,_mm256_set1_ps(1.0/j));
    
    _mm256_store_ps(s,TERM);
    term = (fabs(s[0])>EPS) || (fabs(s[1])>EPS) || (fabs(s[2])>EPS) || (fabs(s[3])>EPS) || 
      (fabs(s[4])>EPS) || (fabs(s[5])>EPS) || (fabs(s[6])>EPS) || (fabs(s[7])>EPS);
    
  }
   
  _mm256_store_ps(s,_mm256_sub_ps(Z,Z0));
  term = (fabs(s[0])<EPS) && (fabs(s[1])<EPS) && (fabs(s[2])<EPS) && (fabs(s[3])<EPS) || 
    (fabs(s[4])<EPS) && (fabs(s[5])<EPS) && (fabs(s[6])<EPS) && (fabs(s[7])<EPS);
  if(term==true)
    return Y;


  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  float alpha;
  bool err = 1;
  Z1 = _mm256_max_ps(Z,C0_BRANCH_PRED);

  __m256 MINUS= _mm256_set1_ps(-1.0);
  __m256 AM2 = ZERO;
  __m256 BM2 = ONE;
  __m256 AM1 = ONE;
  __m256 BM1 = Z1;
  __m256 F   = _mm256_div_ps(AM1,BM1);
  __m256 OLDF= _mm256_set1_ps(INFTY);
  __m256 ALPHA, A, B, INVB, BETA;

  while (err==1 && j<MAXITER+10){
    /* calculate the coefficients of the recursion formulas for j even*/
    alpha = j/2.0;
    ALPHA = _mm256_set1_ps(alpha);
    A = _mm256_add_ps(AM1, _mm256_mul_ps(ALPHA,AM2));
    B = _mm256_add_ps(BM1, _mm256_mul_ps(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_ps(B);
    AM2  = _mm256_mul_ps(AM1,INVB);
    BM2  = _mm256_mul_ps(BM1,INVB);
    AM1  = _mm256_mul_ps(A,INVB);
    BM1  = ONE;

    F    = AM1;
    j++;

    /*calculate the coefficients for j odd*/
    alpha = (j-1.0)/2.0;
    ALPHA = _mm256_set1_ps(alpha);
    BETA  = Z1;
    A = _mm256_add_ps(_mm256_mul_ps(BETA,AM1), _mm256_mul_ps(ALPHA,AM2));
    B = _mm256_add_ps(_mm256_mul_ps(BETA,BM1), _mm256_mul_ps(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_ps(B);
    AM2  = _mm256_mul_ps(AM1,INVB);
    BM2  = _mm256_mul_ps(BM1,INVB);
    AM1  = _mm256_mul_ps(A,INVB);
    BM1  = ONE;

    OLDF = F;
    F    = AM1;
    j++;

    ERR = _mm256_sub_ps(F,OLDF);
    _mm256_store_ps(s,ERR);
    _mm256_store_ps(f,F);
    err = (fabs(s[0]/f[0])>EPS) || (fabs(s[1]/f[1])>EPS) || 
      (fabs(s[2]/f[2])>EPS) || (fabs(s[3]/f[3])>EPS) ||
      (fabs(s[4]/f[4])>EPS) || (fabs(s[5]/f[5])>EPS) || 
      (fabs(s[6]/f[6])>EPS) || (fabs(s[7]/f[7])>EPS);

  }

  /* The output is in F*/
  F = _mm256_mul_ps(F, _mm256_exp_ps(_mm256_mul_ps(MINUS,Z1)));

  /* Return the corresponding answer */
  Y = _mm256_or_ps( _mm256_and_ps(C1,Y),_mm256_andnot_ps(C1,F) );

  return Y;
}


__m256 _mm256_ein_ps(__m256 Z)
{
  int MAXITER = 5;
  float branch_pred     = 3.0162691827353;
  __m256 GAMMA          = _mm256_set1_ps(EULER);
  __m256 C0_BRANCH_PRED = _mm256_set1_ps(branch_pred);
  __m256 ONE            = _mm256_set1_ps(1.0);
  __m256 ZERO           = _mm256_setzero_ps();
  __m256 LOG_P_EULER    = _mm256_add_ps(_mm256_log_ps(Z),GAMMA);
  __m256 ERR,Z0,Z1;
  
  /* branch prediction */
  __m256 C1 = _mm256_cmp_ps(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  float s[8] MEM_ALIGN;
  float f[8] MEM_ALIGN;
  
  
  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = true;
  int j = 1;
  /* If one of the single precision values is large, this loop
   * reaches the MAXITER. Therefore, we set all the values bigger than
   * C0_BRANCH_PRED, to this value. The same is done for the continued 
   * fraction iteration.
   */
  Z0 = _mm256_min_ps(Z,C0_BRANCH_PRED);
  
  __m256 Y = _mm256_sub_ps(ZERO,_mm256_add_ps(GAMMA,_mm256_log_ps(Z0)));
  __m256 PTERM = Z0;
  __m256 TERM  = Z0;
  
  /* The output will be in Y.*/
  while (term==true && j<MAXITER){
    Y     = _mm256_add_ps(Y, TERM);
    j++;
    PTERM = _mm256_mul_ps(Z0,_mm256_mul_ps(PTERM,_mm256_set1_ps(-1.0/j)));
    TERM  = _mm256_mul_ps(PTERM,_mm256_set1_ps(1.0/j));
    
    _mm256_store_ps(s,TERM);
    term = (fabs(s[0])>EPS) || (fabs(s[1])>EPS) || (fabs(s[2])>EPS) || (fabs(s[3])>EPS) || 
      (fabs(s[4])>EPS) || (fabs(s[5])>EPS) || (fabs(s[6])>EPS) || (fabs(s[7])>EPS);
    
  }
  
  /* add the contribution of the log end euler constant */
  Y = _mm256_add_ps( LOG_P_EULER, Y);
  
  _mm256_store_ps(s,_mm256_sub_ps(Z,Z0));
  term = (fabs(s[0])<EPS) && (fabs(s[1])<EPS) && (fabs(s[2])<EPS) && (fabs(s[3])<EPS) &&
    (fabs(s[4])<EPS) && (fabs(s[5])<EPS) && (fabs(s[6])<EPS) && (fabs(s[7])<EPS);
  if(term==true)
    return Y;


  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  float alpha;
  bool err = 1;
  Z1 = _mm256_max_ps(Z,C0_BRANCH_PRED);

  __m256 MINUS= _mm256_set1_ps(-1.0);
  __m256 AM2 = ZERO;
  __m256 BM2 = ONE;
  __m256 AM1 = ONE;
  __m256 BM1 = Z1;
  __m256 F   = _mm256_div_ps(AM1,BM1);
  __m256 OLDF= _mm256_set1_ps(INFTY);
  __m256 ALPHA, A, B, INVB, BETA;

  while (err==1 && j<MAXITER+10){
    /* calculate the coefficients of the recursion formulas for j even*/
    alpha = j/2.0;
    ALPHA = _mm256_set1_ps(alpha);
    A = _mm256_add_ps(AM1, _mm256_mul_ps(ALPHA,AM2));
    B = _mm256_add_ps(BM1, _mm256_mul_ps(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_ps(B);
    AM2  = _mm256_mul_ps(AM1,INVB);
    BM2  = _mm256_mul_ps(BM1,INVB);
    AM1  = _mm256_mul_ps(A,INVB);
    BM1  = ONE;

    F    = AM1;
    j++;

    /*calculate the coefficients for j odd*/
    alpha = (j-1.0)/2.0;
    ALPHA = _mm256_set1_ps(alpha);
    BETA  = Z1;
    A = _mm256_add_ps(_mm256_mul_ps(BETA,AM1), _mm256_mul_ps(ALPHA,AM2));
    B = _mm256_add_ps(_mm256_mul_ps(BETA,BM1), _mm256_mul_ps(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm256_inv_ps(B);
    AM2  = _mm256_mul_ps(AM1,INVB);
    BM2  = _mm256_mul_ps(BM1,INVB);
    AM1  = _mm256_mul_ps(A,INVB);
    BM1  = ONE;

    OLDF = F;
    F    = AM1;
    j++;

    ERR = _mm256_sub_ps(F,OLDF);
    _mm256_store_ps(s,ERR);
    _mm256_store_ps(f,F);
    err = (fabs(s[0]/f[0])>EPS) || (fabs(s[1]/f[1])>EPS) || 
      (fabs(s[2]/f[2])>EPS) || (fabs(s[3]/f[3])>EPS) || 
      (fabs(s[4]/f[4])>EPS) || (fabs(s[5]/f[5])>EPS) || 
      (fabs(s[6]/f[6])>EPS) || (fabs(s[7]/f[7])>EPS);

  }

  F = _mm256_mul_ps(F, _mm256_exp_ps(_mm256_mul_ps(MINUS,Z1)));

  /* The output is in F*/
  /* add the contribution of the log end euler constant */
  F = _mm256_add_ps( LOG_P_EULER, F);

  /* Return the corresponding answer */
  Y = _mm256_or_ps( _mm256_and_ps(C1,Y),_mm256_andnot_ps(C1,F) );

  return Y;
}
#endif //SINGLE


#elif defined __SSE4_2__
#include "math_x86.h"

__m128d _mm_ein_pd(__m128d Z)
{
  int MAXITER = 25;
  __m128d GAMMA          = _mm_set1_pd(EULER);
  __m128d C0_BRANCH_PRED = _mm_set1_pd(3.0162691827353);
  __m128d ONE            = _mm_set1_pd(1.0);
  __m128d ZERO           = _mm_setzero_pd();
  __m128d XEPS           = _mm_set1_pd(EPS);
  __m128d NXEPS          = _mm_set1_pd(-EPS);
  __m128d X100EPS        = _mm_set1_pd(100.0*EPS);
  __m128d NX100EPS       = _mm_set1_pd(-100.0*EPS);
  __m128d LOG_P_EULER    = _mm_add_pd(_mm_log_pd(Z),GAMMA);
  __m128d ERR;

  /* branch prediction */
  __m128d C1 = _mm_cmp_pd(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  double s[2] MEM_ALIGN;

  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = 1;
  int j = 1;
  
  __m128d Y = _mm_sub_pd(ZERO,_mm_add_pd(GAMMA,_mm_log_pd(Z)));
  __m128d PTERM = Z;
  __m128d TERM  = Z;
  
  /* The output will be in Y.*/
  while (term && j<MAXITER){
    Y     = _mm_add_pd(Y, TERM);
    j++;
    PTERM = _mm_mul_pd(Z,_mm_mul_pd(PTERM,_mm_set1_pd(-1.0/j)));
    TERM  = _mm_mul_pd(PTERM,_mm_set1_pd(1.0/j));
    /*since all x are positive*/
    ERR   = _mm_or_pd(_mm_cmp_pd(TERM,XEPS,_CMP_GE_OQ),_mm_cmp_pd(TERM,NXEPS,_CMP_LE_OQ));
    ERR   = _mm_and_pd(ONE, ERR);
    _mm_store_pd(s,ERR);
    term = (sum4(s)>0);
  }

  /* add the contribution of log and Euler constant*/
  Y = _mm_add_pd( LOG_P_EULER, Y);
  if(j<MAXITER)
    return Y;
  
  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  double alpha;
  bool err = 1;
  __m128d MINUS= _mm_set1_pd(-1.0);
  __m128d AM2 = ZERO;
  __m128d BM2 = ONE;
  __m128d AM1 = ONE;
  __m128d BM1 = Z;
  __m128d F   = _mm_div_pd(AM1,BM1);
  __m128d OLDF= _mm_set1_pd(INFTY);
  __m128d ALPHA, A, B, INVB, BETA;

  while (err && j<MAXITER*2){
    /* calculate the coefficients of the recursion formulas for j even*/
    alpha = j/2.0;
    ALPHA = _mm_set1_pd(alpha);
    A = _mm_add_pd(AM1, _mm_mul_pd(ALPHA,AM2));
    B = _mm_add_pd(BM1, _mm_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm_inv_pd(B);
    AM2  = _mm_mul_pd(AM1,INVB);
    BM2  = _mm_mul_pd(BM1,INVB);
    AM1  = _mm_mul_pd(A,INVB);
    BM1  = ONE;

    F    = AM1;
    j++;

    /*calculate the coefficients for j odd*/
    alpha = (j-1.0)/2.0;
    ALPHA = _mm_set1_pd(alpha);
    BETA  = Z;
    A = _mm_add_pd(_mm_mul_pd(BETA,AM1), _mm_mul_pd(ALPHA,AM2));
    B = _mm_add_pd(_mm_mul_pd(BETA,BM1), _mm_mul_pd(ALPHA,BM2));

    /*save new normalized variables for next pass through the loop
     *note: normalization to avoid overflow or underflow*/
    INVB = se_mm_inv_pd(B);
    AM2  = _mm_mul_pd(AM1,INVB);
    BM2  = _mm_mul_pd(BM1,INVB);
    AM1  = _mm_mul_pd(A,INVB);
    BM1  = ONE;

    OLDF = F;
    F    = AM1;
    j++;
    ERR = _mm_mul_pd(_mm_sub_pd(F,OLDF),se_mm_inv_pd(F));
    ERR   = _mm_or_pd(_mm_cmp_pd(ERR,X100EPS,_CMP_GE_OQ),_mm_cmp_pd(ERR,NX100EPS,_CMP_LE_OQ));
    ERR   = _mm_and_pd(ONE, ERR);
    _mm_store_pd(s,ERR);
    err = (sum4(s)>0);
  }

  F = _mm_mul_pd(F, _mm_exp_pd(_mm_mul_pd(MINUS,Z)));
  F = _mm_sub_pd( _mm_add_pd(_mm_log_pd(Z),GAMMA), F);

  /* The output is in F*/
  /* add the contribution of log and Euler constant*/
  F = _mm_add_pd( LOG_P_EULER, F);

  /* Return the corresponding answer */
  Y = _mm_or_pd( _mm_and_pd(C1,Y),_mm_andnot_pd(C1,F) );

  return Y;
}
#endif //AVX

#endif //EXPINT_H_
