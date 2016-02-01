#include "mathint.h"

#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif


double expint(int n, double x)
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


double matlab_expint(double x)
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
	  pterm = -x*pterm/j;
	  term  = pterm/j;
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
	  alpha = j/2.;
	  
	  // calculate A(j), B(j), and f(j)
	  a = am1 + alpha * am2;
	  b = bm1 + alpha * bm2;
   
	  // save new normalized variables for next pass through the loop
	  //  note: normalization to avoid overflow or underflow
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1;
	  
	  f = am1;
	  j++;
   
	  // calculate the coefficients for j odd
	  alpha = (j-1)/2.;
	  beta = x;
	  a = beta  * am1 + alpha * am2;
	  b = beta  * bm1 + alpha * bm2;
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1;
	  oldf= f;
	  f   = am1;
	  j++;
	} // end of while loop
      
      y = exp(-x)*f;

    }


  return y;
}

// expint+log+EULER
double expint_log_euler(double x)
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
//      y = -egamma-log(x);
      y = 0;
      j = 1;
      pterm = x;
      term = x;
      
      while(fabs(term)>EPS)
	{
	  y += term;
	  j ++;
	  pterm = -x*pterm/j;
	  term  = pterm/j;
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
	  alpha = j/2.;
	  
	  // calculate A(j), B(j), and f(j)
	  a = am1 + alpha * am2;
	  b = bm1 + alpha * bm2;
   
	  // save new normalized variables for next pass through the loop
	  //  note: normalization to avoid overflow or underflow
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1;
	  
	  f = am1;
	  j++;
   
	  // calculate the coefficients for j odd
	  alpha = (j-1)/2.;
	  beta = x;
	  a = beta  * am1 + alpha * am2;
	  b = beta  * bm1 + alpha * bm2;
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1;
	  oldf= f;
	  f   = am1;
	  j++;
	} // end of while loop
      
      y = exp(-x)*f;

      egamma=0.57721566490153286061;
      y =  + egamma + log(x) - y;
    }


  return y;
}


// valid for 1<x<=4
static double
cheb_eval_e4(const cheb_series * cs,
	     const double      * X,
	     const double      * q,
	     const int           n
	     )
{
  int j, k, st=0;
  double sum[4] __attribute__((aligned(32)));
  set_zero_4(sum);

#ifdef __AVX__
  __m256d rD1, rDD1, rX1, rY21, rE1, rSS1, rLOG1, rQ1, rS1;
  __m256d rD2, rDD2, rX2, rY22, rE2, rSS2, rLOG2, rQ2, rS2;
  __m256d rD3, rDD3, rX3, rY23, rE3, rSS3, rLOG3, rQ3, rS3;
  __m256d rD4, rDD4, rX4, rY24, rE4, rSS4, rLOG4, rQ4, rS4;

  __m256d ONE = _mm256_set1_pd(1.0);
  __m256d TWO = _mm256_set1_pd(2.0);
  __m256d C1  = _mm256_set1_pd(8.0/3.0);
  __m256d C2  = _mm256_set1_pd(-5.0/3.0);
  __m256d rEULER = _mm256_set1_pd(EULER);

  rS1 = _mm256_setzero_pd();
  rS2 = _mm256_setzero_pd();
  rS3 = _mm256_setzero_pd();
  rS4 = _mm256_setzero_pd();

  for (k=0; k<n/16*16; k+=16)
    {
      rD1  = _mm256_setzero_pd();
      rD2  = _mm256_setzero_pd();
      rD3  = _mm256_setzero_pd();
      rD4  = _mm256_setzero_pd();

      rDD1 = _mm256_setzero_pd();
      rDD2 = _mm256_setzero_pd();
      rDD3 = _mm256_setzero_pd();
      rDD4 = _mm256_setzero_pd();

      rX1  = _mm256_load_pd(X + k   );
      rX2  = _mm256_load_pd(X + k+4 );
      rX3  = _mm256_load_pd(X + k+8 );
      rX4  = _mm256_load_pd(X + k+12);

      // load q
      rQ1  = _mm256_load_pd(q + k   );
      rQ2  = _mm256_load_pd(q + k+4 );
      rQ3  = _mm256_load_pd(q + k+8 );
      rQ4  = _mm256_load_pd(q + k+12);

      // compute log term
      rLOG1 = se_mm256_log_pd(rX1);
      rLOG2 = se_mm256_log_pd(rX2);
      rLOG3 = se_mm256_log_pd(rX3);
      rLOG4 = se_mm256_log_pd(rX4);

      rSS1  = se_mm256_inv_pd(_mm256_mul_pd(rX1,se_mm256_exp_pd(rX1)));
      rSS2  = se_mm256_inv_pd(_mm256_mul_pd(rX2,se_mm256_exp_pd(rX2)));
      rSS3  = se_mm256_inv_pd(_mm256_mul_pd(rX3,se_mm256_exp_pd(rX3)));
      rSS4  = se_mm256_inv_pd(_mm256_mul_pd(rX4,se_mm256_exp_pd(rX4)));

      rX1 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX1)),C2);
      rX2 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX2)),C2);
      rX3 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX3)),C2);
      rX4 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX4)),C2);

      rY21 = _mm256_mul_pd(TWO,rX1);
      rY22 = _mm256_mul_pd(TWO,rX2);
      rY23 = _mm256_mul_pd(TWO,rX3);
      rY24 = _mm256_mul_pd(TWO,rX4);
// FIXME j=24 not 20
      for(j = 24; j>=1; j--) 
	{
	  __m256d rCS;
	  __m256d rTemp1, rTemp2, rTemp3, rTemp4;

	  rCS   = _mm256_set1_pd(cs->c[j]);

	  rTemp1 = _mm256_permute_pd(rD1,5);  // rTemp = rD;
	  rTemp2 = _mm256_permute_pd(rD2,5);  // rTemp = rD;
	  rTemp3 = _mm256_permute_pd(rD3,5);  // rTemp = rD;
	  rTemp4 = _mm256_permute_pd(rD4,5);  // rTemp = rD;

	  rD1   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY21,rD1),rDD1));
	  rD2   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY22,rD2),rDD2));
	  rD3   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY23,rD3),rDD3));
	  rD4   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY24,rD4),rDD4));

	  rDD1 = _mm256_permute_pd(rTemp1,5);  // rDD = rTemp;
	  rDD2 = _mm256_permute_pd(rTemp2,5);  // rDD = rTemp;
	  rDD3 = _mm256_permute_pd(rTemp3,5);  // rDD = rTemp;
	  rDD4 = _mm256_permute_pd(rTemp4,5);  // rDD = rTemp;
	}
      
      {
	__m256d rCS;
	__m256d rTemp1, rTemp2, rTemp3, rTemp4;

	rCS   = _mm256_set1_pd(0.5*cs->c[0]);

	rTemp1 = _mm256_permute_pd(rD1,5);  // rTemp = rD;
	rTemp2 = _mm256_permute_pd(rD2,5);  // rTemp = rD;
	rTemp3 = _mm256_permute_pd(rD3,5);  // rTemp = rD;
	rTemp4 = _mm256_permute_pd(rD4,5);  // rTemp = rD;
	rD1    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX1,rD1),rDD1),rCS);
	rD2    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX2,rD2),rDD2),rCS);
	rD3    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX3,rD3),rDD3),rCS);
	rD4    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX4,rD4),rDD4),rCS);
      }

      rE1 = _mm256_mul_pd(rSS1,_mm256_add_pd(ONE,rD1));
      rE2 = _mm256_mul_pd(rSS2,_mm256_add_pd(ONE,rD2));
      rE3 = _mm256_mul_pd(rSS3,_mm256_add_pd(ONE,rD3));
      rE4 = _mm256_mul_pd(rSS4,_mm256_add_pd(ONE,rD4));

      // add log term and EULER constant
      rS1 = _mm256_add_pd(rS1,_mm256_mul_pd(rQ1,_mm256_add_pd(rE1,_mm256_add_pd(rLOG1,rEULER))));
      rS2 = _mm256_add_pd(rS2,_mm256_mul_pd(rQ2,_mm256_add_pd(rE2,_mm256_add_pd(rLOG2,rEULER))));
      rS3 = _mm256_add_pd(rS3,_mm256_mul_pd(rQ3,_mm256_add_pd(rE3,_mm256_add_pd(rLOG3,rEULER))));
      rS4 = _mm256_add_pd(rS4,_mm256_mul_pd(rQ4,_mm256_add_pd(rE4,_mm256_add_pd(rLOG4,rEULER))));

    }
  rS1 = _mm256_add_pd(rS1,_mm256_add_pd(rS2,_mm256_add_pd(rS3,rS4)));
  _mm256_store_pd(sum, rS1);
  sum[0] += sum[1]+sum[2]+sum[3];
  st = k;
#endif
  double s,x;
  for (k=st; k<n; k++)
    {
      double d  = 0.0;
      double dd = 0.0;
      x = X[k];
      s = 1.0/x*exp(-x);
      x = (8.0/x-5.0)/3.0;

      double y2 = 2.0 * x;
      // FIXME 24 not 20
      for(j = 24; j>=1; j--) 
	{
	  double temp = d;
	  d = y2*d - dd + cs->c[j];
	  dd = temp;
	}
      
      {
	double temp = d;
	d = x*d - dd + 0.5 * cs->c[0];
      }

      // add log term and EULER constant
      sum[0] += q[k]*(s*(1.0+d)+log(X[k])+EULER);
    }
  return sum[0];
}

// valid for 4<x<=32
static double
cheb_eval_e32(const cheb_series * cs,
	      const double      * X,
	      const double      * q,
	      const int           n
	      )
{
  int j, k, st=0;
  double sum[4] __attribute__((aligned(32)));
  set_zero_4(sum);

#ifdef __AVX__
  __m256d rD1, rDD1, rX1, rY21, rE1, rSS1, rLOG1, rQ1, rS1;
  __m256d rD2, rDD2, rX2, rY22, rE2, rSS2, rLOG2, rQ2, rS2;
  __m256d rD3, rDD3, rX3, rY23, rE3, rSS3, rLOG3, rQ3, rS3;
  __m256d rD4, rDD4, rX4, rY24, rE4, rSS4, rLOG4, rQ4, rS4;

  __m256d ONE    = _mm256_set1_pd(1.0);
  __m256d TWO    = _mm256_set1_pd(2.0);
  __m256d C1     = _mm256_set1_pd(8.0);
  __m256d C2     = _mm256_set1_pd(-1.0);
  __m256d rEULER = _mm256_set1_pd(EULER);

  rS1 = _mm256_setzero_pd();
  rS2 = _mm256_setzero_pd();
  rS3 = _mm256_setzero_pd();
  rS4 = _mm256_setzero_pd();

  for (k=0; k<n/16*16; k+=16)
    {
      rD1  = _mm256_setzero_pd();
      rD2  = _mm256_setzero_pd();
      rD3  = _mm256_setzero_pd();
      rD4  = _mm256_setzero_pd();

      rDD1 = _mm256_setzero_pd();
      rDD2 = _mm256_setzero_pd();
      rDD3 = _mm256_setzero_pd();
      rDD4 = _mm256_setzero_pd();

      rX1  = _mm256_load_pd(X + k   );
      rX2  = _mm256_load_pd(X + k+4 );
      rX3  = _mm256_load_pd(X + k+8 );
      rX4  = _mm256_load_pd(X + k+12);

      // load q
      rQ1  = _mm256_load_pd(q + k   );
      rQ2  = _mm256_load_pd(q + k+4 );
      rQ3  = _mm256_load_pd(q + k+8 );
      rQ4  = _mm256_load_pd(q + k+12);

      // compute log term
      rLOG1 = se_mm256_log_pd(rX1);
      rLOG2 = se_mm256_log_pd(rX2);
      rLOG3 = se_mm256_log_pd(rX3);
      rLOG4 = se_mm256_log_pd(rX4);

      rSS1  = se_mm256_inv_pd(_mm256_mul_pd(rX1,se_mm256_exp_pd(rX1)));
      rSS2  = se_mm256_inv_pd(_mm256_mul_pd(rX2,se_mm256_exp_pd(rX2)));
      rSS3  = se_mm256_inv_pd(_mm256_mul_pd(rX3,se_mm256_exp_pd(rX3)));
      rSS4  = se_mm256_inv_pd(_mm256_mul_pd(rX4,se_mm256_exp_pd(rX4)));

      rX1 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX1)),C2);
      rX2 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX2)),C2);
      rX3 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX3)),C2);
      rX4 = _mm256_add_pd(_mm256_mul_pd(C1,se_mm256_inv_pd(rX4)),C2);

      rY21 = _mm256_mul_pd(TWO,rX1);
      rY22 = _mm256_mul_pd(TWO,rX2);
      rY23 = _mm256_mul_pd(TWO,rX3);
      rY24 = _mm256_mul_pd(TWO,rX4);

      for(j = 25; j>=1; j--) 
	{
	  __m256d rCS;
	  __m256d rTemp1, rTemp2, rTemp3, rTemp4;

	  rCS   = _mm256_set1_pd(cs->c[j]);

	  rTemp1 = _mm256_permute_pd(rD1,5);  // rTemp = rD;
	  rTemp2 = _mm256_permute_pd(rD2,5);  // rTemp = rD;
	  rTemp3 = _mm256_permute_pd(rD3,5);  // rTemp = rD;
	  rTemp4 = _mm256_permute_pd(rD4,5);  // rTemp = rD;

	  rD1   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY21,rD1),rDD1));
	  rD2   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY22,rD2),rDD2));
	  rD3   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY23,rD3),rDD3));
	  rD4   = _mm256_add_pd(rCS,_mm256_sub_pd(_mm256_mul_pd(rY24,rD4),rDD4));

	  rDD1 = _mm256_permute_pd(rTemp1,5);  // rDD = rTemp;
	  rDD2 = _mm256_permute_pd(rTemp2,5);  // rDD = rTemp;
	  rDD3 = _mm256_permute_pd(rTemp3,5);  // rDD = rTemp;
	  rDD4 = _mm256_permute_pd(rTemp4,5);  // rDD = rTemp;
	}
      
      {
	__m256d rCS;
	__m256d rTemp1, rTemp2, rTemp3, rTemp4;

	rCS   = _mm256_set1_pd(0.5*cs->c[0]);

	rTemp1 = _mm256_permute_pd(rD1,5);  // rTemp = rD;
	rTemp2 = _mm256_permute_pd(rD3,5);  // rTemp = rD;
	rTemp3 = _mm256_permute_pd(rD3,5);  // rTemp = rD;
	rTemp4 = _mm256_permute_pd(rD4,5);  // rTemp = rD;

	rD1    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX1,rD1),rDD1),rCS);
	rD2    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX2,rD2),rDD2),rCS);
	rD3    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX3,rD3),rDD3),rCS);
	rD4    = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(rX4,rD4),rDD4),rCS);
      }

      rE1 = _mm256_mul_pd(rSS1,_mm256_add_pd(ONE,rD1));
      rE2 = _mm256_mul_pd(rSS2,_mm256_add_pd(ONE,rD2));
      rE3 = _mm256_mul_pd(rSS3,_mm256_add_pd(ONE,rD3));
      rE4 = _mm256_mul_pd(rSS4,_mm256_add_pd(ONE,rD4));

      // add log term and EULER constant
      rS1 = _mm256_add_pd(rS1,_mm256_mul_pd(rQ1,_mm256_add_pd(rE1,_mm256_add_pd(rLOG1,rEULER))));
      rS2 = _mm256_add_pd(rS2,_mm256_mul_pd(rQ2,_mm256_add_pd(rE2,_mm256_add_pd(rLOG2,rEULER))));
      rS3 = _mm256_add_pd(rS3,_mm256_mul_pd(rQ3,_mm256_add_pd(rE3,_mm256_add_pd(rLOG3,rEULER))));
      rS4 = _mm256_add_pd(rS4,_mm256_mul_pd(rQ4,_mm256_add_pd(rE4,_mm256_add_pd(rLOG4,rEULER))));
      
    }
  rS1 = _mm256_add_pd(rS1,_mm256_add_pd(rS2,_mm256_add_pd(rS3,rS4)));
  _mm256_store_pd(sum, rS1);

  sum[0] += sum[1]+sum[2]+sum[3];
  st = k;
#endif
  double s,x;
  for (k=st; k<n; k++)
    {
      double d  = 0.0;
      double dd = 0.0;
      x = X[k];
      s = 1.0/x*exp(-x);
      x = 8.0/x-1.0;
	
      double y2 = 2.0 * x;
      for(j = 25; j>=1; j--) 
	{
	  double temp = d;
	  d = y2*d - dd + cs->c[j];
	  dd = temp;
	}
	
      {
	double temp = d;
	d = x*d - dd + 0.5 * cs->c[0];
      }
      // add log term and EULER constant
      sum[0] += q[k]*(s*(1.0+d)+log(X[k])+EULER);
    }
  return sum[0];
}

static double
cheb_eval_eMAX(const double   * X, 
	       const double   * q,
	       const int        n
	       )
{
  int k,st=0;
  double sum[4] __attribute__((aligned(32)));
  set_zero_4(sum);
#ifdef __AVX__
  __m256d rS1,rS2, rS3, rS4;
  __m256d rX1, rX2, rX3, rX4, rEULER;
  __m256d rQ1,rQ2, rQ3, rQ4;
  rEULER = _mm256_set1_pd(EULER);
  rS1 = _mm256_setzero_pd();
  rS2 = _mm256_setzero_pd();
  rS3 = _mm256_setzero_pd();
  rS4 = _mm256_setzero_pd();

  for (k=0; k<n/16*16; k+=16)
    {  
      // load X 
      rX1 = _mm256_load_pd( X + k );
      rX2 = _mm256_load_pd( X + k+4 );
      rX3 = _mm256_load_pd( X + k+8 );
      rX4 = _mm256_load_pd( X + k+12 );
      
      // load q
      rQ1 = _mm256_load_pd( q + k );
      rQ2 = _mm256_load_pd( q + k+4 );
      rQ3 = _mm256_load_pd( q + k+8 );
      rQ4 = _mm256_load_pd( q + k+12 );

      // add log and Euler constants to zero as expint=0 for x>32
      rS1 = _mm256_add_pd(rS1,_mm256_mul_pd(rQ1,_mm256_add_pd(se_mm256_log_pd(rX1),rEULER)));
      rS2 = _mm256_add_pd(rS2,_mm256_mul_pd(rQ2,_mm256_add_pd(se_mm256_log_pd(rX2),rEULER)));
      rS3 = _mm256_add_pd(rS3,_mm256_mul_pd(rQ3,_mm256_add_pd(se_mm256_log_pd(rX3),rEULER)));
      rS4 = _mm256_add_pd(rS4,_mm256_mul_pd(rQ4,_mm256_add_pd(se_mm256_log_pd(rX4),rEULER)));

    }

  rS1 = _mm256_add_pd(rS1,_mm256_add_pd(rS2,_mm256_add_pd(rS3,rS4)));
  _mm256_store_pd(sum, rS1);
  sum[0] += sum[1]+sum[2]+sum[3];
  st = k;
#endif
  for (k=st; k<n; k++)
    sum[0] += q[k] *(log(X[k])+EULER);
  return sum[0];
}


static double E12_data[16] = {
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
static cheb_series E12_cs = {
  E12_data,
  15
};

static double AE13_data[25] = {
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
static cheb_series AE13_cs = {
  AE13_data,
  24
};

static double AE14_data[26] = {
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
static cheb_series AE14_cs = {
  AE14_data,
  25
};

double gsl_expint_log_euler_1(const double *x, const double * q, const int n)
{
  /* CHECK_POINTER(result) */
  if(n==0) return 0;
  return cheb_eval_e1(&E12_cs, x, q, n);
}

double gsl_expint_log_euler_4(const double *x, const double * q, const int n)
{
  /* CHECK_POINTER(result) */
  if(n==0) return 0;
  return cheb_eval_e4(&AE13_cs, x, q,n);
}

double gsl_expint_log_euler_32(const double *x, const double * q, const int n)
{
  /* CHECK_POINTER(result) */
  if(n==0) return 0;
  return cheb_eval_e32(&AE14_cs, x, q,n);
}

double gsl_expint_log_euler_MAX(const double *x, const double * q, const int n)
{
  /* CHECK_POINTER(result) */
  if(n==0) return 0;
  return cheb_eval_eMAX(x, q, n);
}

void
constants(double *t, double *E, double *T, double *W, int M, const double loga, const double logb, const int gl_order)
{
  // non uniform tabulated expints
  double c;
  double h = (logb-loga)/(M-1.);

#ifdef _OPENMP
#pragma omp parallel for shared (t,E)
#endif
  for (int k=0; k<M; k++)
    {
      c = loga + k*h;
      t[k] = pow(10.,c);
      if(t[k]>66)
	E[k] = 0;
      else
	E[k] = matlab_expint(t[k]);
    }
 
  // tabulated G-L quadrature nodes and weights
  if(gl_order == 8)
    {
      T[0] = 0.960289856497536398194370121928;
      T[1] = 0.796666477413626727965834106726;
      T[2] = 0.525532409916328990817646626965;
      T[3] = 0.183434642495649835591819964975;
      T[4] = -0.183434642495649835591819964975;
      T[5] = -0.525532409916328990817646626965;
      T[6] = -0.796666477413626727965834106726;
      T[7] = -0.960289856497536398194370121928;
      
      W[0] = 0.101228536290376813777669440242;
      W[1] = 0.222381034453374426540506192396;
      W[2] = 0.313706645877887324580512995453;
      W[3] = 0.362683783378362101235126147003;
      W[4] = 0.362683783378362101235126147003;
      W[5] = 0.313706645877887324580512995453;
      W[6] = 0.222381034453374426540506192396;
      W[7] = 0.101228536290376813777669440242;
    }
  else if(gl_order==4)
    {
      T[0] = 0.861136311594052017426292877644;
      T[1] = 0.339981043584855979755587895852;
      T[2] = -0.339981043584855979755587895852;
      T[3] = -0.861136311594052017426292877644;
      W[0] = 0.347854845137454016246181254246;
      W[1] = 0.652145154862545983753818745754;
      W[2] = 0.652145154862545983753818745754;
      W[3] = 0.347854845137454016246181254246;
    }
  else if (gl_order == 2)
    {
      T[0] = -0.577350269189625731058868041146;
      T[1] = 0.577350269189625731058868041146;
      W[0] = 1.0;
      W[1] = 1.0;
    }
  else
    glwt_fast(gl_order, -1., 1., T,W);

}

inline void
mapping(double a, double b, double *T1, double *W1, double* T, double *W, const int gl_order)
{

  for (int k=0; k<gl_order; k++)
    {
      T[k] = (a+b)/2. + (b-a)/2.*T1[k];
      W[k] = (b-a)/2.*W1[k];
    }
}

inline void
find_idx(double x, int *idx, double *z, double h, const double loga, const double logb)
{
  // we assume that the interval is [10^loga, 10^logb]
  *z = log10(x);
  *idx = ceil((*z-loga)/h);
}

static inline 
double glquad(double a, double b, double* T, double* W, int gl_order)
{
  double s=0;
  for (int k=0; k<gl_order; k++)
    {
      s += W[k] * exp(-T[k])/T[k];
    }
  return s;
}

double 
k0_term(double * t,  double * E,
	double * T1, double * W1,
	const double * restrict x, const double * restrict q,
	int      N,  int      M,
	const int gl_order,
	const double loga, const double logb
	)
{
  int i, st=0;
  double s=0.;
  
#ifdef __AVX__
  double idx[8]  __attribute__((aligned(32)));
  double sum[4]  __attribute__((aligned(32)));
  const __m256d ONE_LOG10 = _mm256_set1_pd(0.434294481903251761156781185491);
  const __m256d LOGA      = _mm256_set1_pd(-loga);
  const __m256d ONE_H     = _mm256_set1_pd((M-1.)/(logb-loga));
  const __m256d HALF      = _mm256_set1_pd(0.5);
  const __m256d rEULER    = _mm256_set1_pd(EULER);


  __m256d rX1, rZ1, rLOGX1, rB1, IDX1, rINT1, rS;
  __m256d rX2, rZ2, rLOGX2, rB2, IDX2, rINT2;
  __m256d rBPA1, rBMA1, rW1, rT1, rQ1;
  __m256d rBPA2, rBMA2, rW2, rT2, rQ2;
  rS  = _mm256_setzero_pd();
  for (i=0; i<N/8*8; i+=8)
    {
      // find the index idx and load corresponding data
      rX1    = _mm256_load_pd(x+i);
      rX2    = _mm256_load_pd(x+i+4);
      rLOGX1 = se_mm256_log_pd(rX1);
      rLOGX2 = se_mm256_log_pd(rX2);
      rZ1    = _mm256_mul_pd(rLOGX1,ONE_LOG10);  // log10(x) = log(x)* (1/log(10))
      rZ2    = _mm256_mul_pd(rLOGX2,ONE_LOG10);  // log10(x) = log(x)* (1/log(10))
      IDX1   = _mm256_ceil_pd(_mm256_mul_pd(_mm256_add_pd(rZ1,LOGA),ONE_H));
      IDX2   = _mm256_ceil_pd(_mm256_mul_pd(_mm256_add_pd(rZ2,LOGA),ONE_H));
      
      _mm256_store_pd(idx,IDX1);
      _mm256_store_pd(idx+4,IDX2);
      
      rB1    = _mm256_set_pd(t[(int)idx[3]], t[(int)idx[2]], t[(int)idx[1]], t[(int)idx[0]]);
      rB2    = _mm256_set_pd(t[(int)idx[7]], t[(int)idx[6]], t[(int)idx[5]], t[(int)idx[4]]);

      rBPA1  = _mm256_mul_pd(HALF,_mm256_add_pd(rB1, rX1));
      rBPA2  = _mm256_mul_pd(HALF,_mm256_add_pd(rB2, rX2));
      rBMA1  = _mm256_mul_pd(HALF,_mm256_sub_pd(rB1, rX1));
      rBMA2  = _mm256_mul_pd(HALF,_mm256_sub_pd(rB2, rX2));

      rINT1  = _mm256_setzero_pd();
      rINT2  = _mm256_setzero_pd();


      // mapping and integration
      // FIXME: gl_order should be even.
      // FIXME: use for loop if gl_order is greater than 2 and remove "int k=0"
      for (int k=0; k<gl_order; k+=2)
      //int k=0;
	{
	  // mapping
	  __m256d rTMP, rF1, rF2;
	  rTMP  = _mm256_set1_pd(W1[k]);
	  rW1   = _mm256_mul_pd(rTMP,rBMA1);
	  rW2   = _mm256_mul_pd(rTMP,rBMA2);
	  rTMP  = _mm256_set1_pd(T1[k]);
	  rT1   = _mm256_add_pd(rBPA1, _mm256_mul_pd(rBMA1,rTMP));
	  rT2   = _mm256_add_pd(rBPA2, _mm256_mul_pd(rBMA2,rTMP));

	  // integration 
	  rF1 = se_mm256_inv_pd(_mm256_mul_pd(se_mm256_exp_pd(rT1),rT1));
	  rF2 = se_mm256_inv_pd(_mm256_mul_pd(se_mm256_exp_pd(rT2),rT2));
	  rINT1 = _mm256_add_pd(rINT1, _mm256_mul_pd(rW1,rF1));
	  rINT2 = _mm256_add_pd(rINT2, _mm256_mul_pd(rW2,rF2));

	  rTMP  = _mm256_set1_pd(W1[k+1]);
	  rW1   = _mm256_mul_pd(rTMP,rBMA1);
	  rW2   = _mm256_mul_pd(rTMP,rBMA2);
	  rTMP  = _mm256_set1_pd(T1[k+1]);
	  rT1   = _mm256_add_pd(rBPA1, _mm256_mul_pd(rBMA1,rTMP));
	  rT2   = _mm256_add_pd(rBPA2, _mm256_mul_pd(rBMA2,rTMP));

	  // integration 
	  rF1 = se_mm256_inv_pd(_mm256_mul_pd(se_mm256_exp_pd(rT1),rT1));
	  rF2 = se_mm256_inv_pd(_mm256_mul_pd(se_mm256_exp_pd(rT2),rT2));
	  rINT1 = _mm256_add_pd(rINT1, _mm256_mul_pd(rW1,rF1));
	  rINT2 = _mm256_add_pd(rINT2, _mm256_mul_pd(rW2,rF2));

	}

      rB1 = _mm256_set_pd(E[(int)idx[3]], E[(int)idx[2]], E[(int)idx[1]], E[(int)idx[0]]);
      rB2 = _mm256_set_pd(E[(int)idx[7]], E[(int)idx[6]], E[(int)idx[5]], E[(int)idx[4]]);
      rB1 = _mm256_add_pd(rINT1, _mm256_add_pd(rB1,_mm256_add_pd(rEULER,rLOGX1)));
      rB2 = _mm256_add_pd(rINT2, _mm256_add_pd(rB2,_mm256_add_pd(rEULER,rLOGX2)));

      rQ1 = _mm256_load_pd(q+i);
      rQ2 = _mm256_load_pd(q+i+4);

      rS = _mm256_add_pd(rS,_mm256_mul_pd(rQ1,rB1));
      rS = _mm256_add_pd(rS,_mm256_mul_pd(rQ2,rB2));
    }
  _mm256_store_pd(sum,rS);
  s += sum[0] + sum[1] + sum[2] + sum[3];
  
  st = i; // to do the reminder of the particles
#endif
  double *T = (double*) _mm_malloc(gl_order*sizeof(double),32);
  double *W = (double*) _mm_malloc(gl_order*sizeof(double),32);

  int iidx;
  double z,a,b;
  for (i=st; i<N; i++)
    {
      find_idx(x[i], &iidx, &z, (logb-loga)/(M-1.), loga, logb);
      a = x[i]; b = t[iidx];
      mapping(a, b, T1, W1, T, W, gl_order);
      if(x[i]>32)
	s += q[i]*(log(x[i])+EULER);
      else
	s += q[i]*(glquad(a, b, T, W, gl_order) + E[iidx]+log(x[i])+EULER);

    }
  _mm_free(T); _mm_free(W);

  return s;

}

