#ifndef EXPINT_H
#define EXPINT_H
/*------------------------------------------------------------------------
 *This function computes expint, based on the MATLAB expint implementation
 *------------------------------------------------------------------------
 */
#include "stdbool.h"

#define EPS    1.110223024625157e-16
#define EULER  0.577215664901532860606512090082402431042
#define INFTY  1.0e308

/* sum up 4 elements in a vector, used for SIMD op.*/
double sum4(double*p){
  /*Unexpected behaviour if p has more or less
   *elements than 4 */
  return (p[0]+p[1]+p[2]+p[3]);
}

double polyval(double *p, double x)
{
  int i;
  double y=p[0];
  for (i=1;i<9;i++)
    {
      y = y*x + p[i];
    }
  return y;
}

double expint(double x)
{
  int j;
  double am1,am2,bm1,bm2,f,oldf,pterm,term,y,polyv,egamma,a,b,alpha,beta;
  double INF = 1.e308;
  double p[] = {-3.602693626336023e-09, -4.819538452140960e-07, -2.569498322115933e-05,
		-6.973790859534190e-04, -1.019573529845792e-02, -7.811863559248197e-02,
		-3.012432892762715e-01, -7.773807325735529e-01,  8.267661952366478e+00};

  polyv = polyval(p,x);
  if(polyv>=0)
    {
      egamma=0.57721566490153286061;
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

/*------------------------------------------------------------------------
 *This function computes expint+log+euler, based on MATLAB expint implementation
 *------------------------------------------------------------------------
 */
double expint_log_euler(double x)
{
  double egamma = 0.57721566490153286061;
  if(fabs(x)<1.0e-16)
    printf("WTF\n");
  if(fabs(x)<1.0e-16)
    return 0.0;
  else
    return expint(x)+log(x)+egamma;
}


double expints(int n, double x)
{
  int i,ii,nm1, MAXIT = 300;
  double a,b,c,d,del,fact,h,psi,ans=0;
  double FPMIN = 1.e-308;

  nm1=n-1;
  if (n < 0 || x < 0.0 || (fabs(x)<1e-18 && (n==0 || n==1)))
    {printf("bad arguments in expint %f\n",x);return -1;}
  else {
    if (n == 0) ans=exp(-x)/x;
    else {
      if (x == 0.0) ans=1.0/nm1;

      else {
	if (x > 1.0) {
	  b=x+n;
	  c=1.0/FPMIN;
	  d=1.0/b;
	  h=d;
	  for (i=1;i<=MAXIT;i++) {
	    a = -i*(nm1+i);
	    b += 2.0;
	    d=1.0/(a*d+b);
	    c=b+a/c;
	    del=c*d;
	    h *= del;
	    if (fabs(del-1.0) < 2.*EPS) {
	      ans=h*exp(-x);
	      return ans;
	    }
	  }
	  printf("continued fraction failed in expint");
	} else {
	  ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
	  fact=1.0;
	  for (i=1;i<=MAXIT;i++) {
	    fact *= -x/i;
	    if (i != nm1) del = -fact/(i-nm1);
	    else {
	      psi = -EULER;
	      for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
	      del=fact*(-log(x)+psi);
	    }
	    ans += del;
	    if (fabs(del) < fabs(ans)*EPS) return ans;
	  }
	  printf("series failed in expint");
	}
      }
    }
  }
  return ans;
}

#ifdef __AVX__
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
__m256d _mm256_ein_pd(__m256d Z)
{
  int MAXITER = 25;
  __m256d GAMMA          = _mm256_set1_pd(0.57721566490153286061);
  __m256d C0_BRANCH_PRED = _mm256_set1_pd(3.0162691827353);
  __m256d XEPS           = _mm256_set1_pd(EPS);
  __m256d NXEPS          = _mm256_set1_pd(-EPS);
  __m256d X100EPS        = _mm256_set1_pd(100.0*EPS);
  __m256d NX100EPS       = _mm256_set1_pd(-100.0*EPS);
  __m256d ONE            = _mm256_set1_pd(1.0);
  __m256d ZERO           = _mm256_setzero_pd();
  __m256d LOG_P_EULER    = _mm256_add_pd(se_mm256_log_pd(Z),GAMMA);
  __m256d ERR;
  
  /* branch prediction */
  __m256d C1 = _mm256_cmp_pd(Z,C0_BRANCH_PRED,_CMP_LT_OQ);
  
  double s[4] MEM_ALIGN;
  
  /*If polyv >=0, Z < 3.0162691827353.*/
  bool term = true;
  int j = 1;
  
  __m256d Y = _mm256_sub_pd(ZERO,_mm256_add_pd(GAMMA,se_mm256_log_pd(Z)));
  __m256d PTERM = Z;
  __m256d TERM  = Z;
  
  /* The output will be in Y.*/
  while (term && j<MAXITER){
    Y     = _mm256_add_pd(Y, TERM);
    j++;
    PTERM = _mm256_mul_pd(Z,_mm256_mul_pd(PTERM,_mm256_set1_pd(-1.0/j)));
    TERM  = _mm256_mul_pd(PTERM,_mm256_set1_pd(1.0/j));

    ERR   = _mm256_or_pd(_mm256_cmp_pd(TERM,XEPS,_CMP_GE_OQ),_mm256_cmp_pd(TERM,NXEPS,_CMP_LE_OQ));
    ERR   = _mm256_and_pd(ONE, ERR);
    _mm256_store_pd(s,ERR);
    term = (sum4(s)>0);
  }
  /* add the contribution of the log end euler constant */
  Y = _mm256_add_pd( LOG_P_EULER, Y);

  if(j<MAXITER)
    return Y;

  /* If polyv <0, Z >= 3.0162691827353. 
   *This is the continued fraction.  */
  j = 2;
  double alpha;
  bool err = 1;
  __m256d MINUS= _mm256_set1_pd(-1.0);
  __m256d AM2 = ZERO;
  __m256d BM2 = ONE;
  __m256d AM1 = ONE;
  __m256d BM1 = Z;
  __m256d F   = _mm256_div_pd(AM1,BM1);
  __m256d OLDF= _mm256_set1_pd(INFTY);
  __m256d ALPHA, A, B, INVB, BETA;

  while (err && j<MAXITER*2){
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
    BETA  = Z;
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
    ERR = _mm256_mul_pd(_mm256_sub_pd(F,OLDF),se_mm256_inv_pd(F));
    ERR   = _mm256_or_pd(_mm256_cmp_pd(ERR,X100EPS,_CMP_GE_OQ),_mm256_cmp_pd(ERR,NX100EPS,_CMP_LE_OQ));
    ERR   = _mm256_and_pd(ONE, ERR);

    _mm256_store_pd(s,ERR);
    err = (sum4(s)>0);
  }

  F = _mm256_mul_pd(F, se_mm256_exp_pd(_mm256_mul_pd(MINUS,Z)));

  /* The output is in F*/
  /* Return the corresponding answer */
  F = _mm256_add_pd( LOG_P_EULER, F);
  Y = _mm256_or_pd( _mm256_and_pd(C1,Y),_mm256_andnot_pd(C1,F) );

  return Y;
}

#elif defined __SSE4_2__
#include "math_x86.h"

__m128d _mm_ein_pd(__m128d Z)
{
  int MAXITER = 25;
  __m128d GAMMA          = _mm_set1_pd(0.57721566490153286061);
  __m128d C0_BRANCH_PRED = _mm_set1_pd(3.0162691827353);
  __m128d ONE            = _mm_set1_pd(1.0);
  __m128d ZERO           = _mm_setzero_pd();
  __m128d XEPS           = _mm_set1_pd(EPS);
  __m128d NXEPS          = _mm_set1_pd(-EPS);
  __m128d X100EPS        = _mm_set1_pd(100.0*EPS);
  __m128d NX100EPS       = _mm_set1_pd(-100.0*EPS);
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
  /* Return the corresponding answer */
  Y = _mm_or_pd( _mm_and_pd(C1,Y),_mm_andnot_pd(C1,F) );
  Y = _mm_add_pd( _mm_add_pd(_mm_log_pd(Z),GAMMA), Y);
  return Y;
}
#endif //AVX
#endif //EXPINT_H_
