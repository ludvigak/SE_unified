#ifndef EXPINT_H
#define EXPINT_H
/*------------------------------------------------------------------------
 *This function computes expint, based on the MATLAB expint implementation
 *------------------------------------------------------------------------
 */
#define EPS    1.110223024625157e-16
#define EULER  0.577215664901532860606512090082402431042
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
    return expint(x)+log(x)+EULER;
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
	  
	  /* calculate A(j), B(j), and f(j)*/
	  a = am1 + alpha * am2;
	  b = bm1 + alpha * bm2;
   
	  /* save new normalized variables for next pass through the loop
	   * note: normalization to avoid overflow or underflow*/
	  am2 = am1 / b;
	  bm2 = bm1 / b;
	  am1 = a   / b;
	  bm1 = 1;
	  
	  f = am1;
	  j++;
   
	  /* calculate the coefficients for j odd*/
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
	} /* end of while loop*/
      
      y = exp(-x)*f;

      egamma=0.57721566490153286061;
      y =  + egamma + log(x) - y;
    }

  return y;
}
#endif
