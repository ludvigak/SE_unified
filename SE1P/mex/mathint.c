#include "mathint.h"

#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif


double GL_value[] = {
  0.0000875560264341,
  0.0004612700113121,
  0.0011333756872430,
  0.0021036207325094,
  0.0033714435498935,
  0.0049360907541328,
  0.0067966286377069,
  0.0089519457821408,
  0.0114007542680463,
  0.0141415906264317,
  0.0171728167840174,
  0.0204926210731500,
  0.0240990193293678,
  0.0279898560848899,
  0.0321628058610418,
  0.0366153745605260,
  0.0413449009595198,
  0.0463485582991216,
  0.0516233559754209,
  0.0571661413273014,
  0.0629736015209841,
  0.0690422655302257,
  0.0753685062110155,
  0.0819485424695466,
  0.0887784415221781,
  0.0958541212460431,
  0.1031713526189034,
  0.1107257622467940,
  0.1185128349779526,
  0.1265279166014690,
  0.1347662166290456,
  0.1432228111582063,
  0.1518926458152429,
  0.1607705387761404,
  0.1698511838636770,
  0.1791291537188462,
  0.1885989030447076,
  0.1982547719207257,
  0.2080909891856185,
  0.2181016758866909,
  0.2282808487935948,
  0.2386224239744122,
  0.2491202204319278,
  0.2597679637979140,
  0.2705592900832239,
  0.2814877494814479,
  0.2925468102238625,
  0.3037298624833663,
  0.3150302223250705,
  0.3264411357011823,
  0.3379557824877933,
  0.3495672805611614,
  0.3612686899110478,
  0.3730530167886529,
  0.3849132178866700,
  0.3968422045489604,
  0.4088328470073314,
  0.4208779786428876,
  0.4329704002694061,
  0.4451028844361781,
  0.4572681797477423,
  0.4694590151979302,
  0.4816681045156332,
  0.4938881505196921,
  0.5061118494803079,
  0.5183318954843668,
  0.5305409848020698,
  0.5427318202522577,
  0.5548971155638218,
  0.5670295997305939,
  0.5791220213571124,
  0.5911671529926685,
  0.6031577954510396,
  0.6150867821133300,
  0.6269469832113471,
  0.6387313100889522,
  0.6504327194388386,
  0.6620442175122067,
  0.6735588642988177,
  0.6849697776749295,
  0.6962701375166337,
  0.7074531897761375,
  0.7185122505185521,
  0.7294407099167761,
  0.7402320362020860,
  0.7508797795680722,
  0.7613775760255878,
  0.7717191512064052,
  0.7818983241133091,
  0.7919090108143816,
  0.8017452280792743,
  0.8114010969552925,
  0.8208708462811538,
  0.8301488161363231,
  0.8392294612238596,
  0.8481073541847571,
  0.8567771888417937,
  0.8652337833709545,
  0.8734720833985310,
  0.8814871650220474,
  0.8892742377532059,
  0.8968286473810967,
  0.9041458787539569,
  0.9112215584778219,
  0.9180514575304535,
  0.9246314937889845,
  0.9309577344697743,
  0.9370263984790159,
  0.9428338586726985,
  0.9483766440245791,
  0.9536514417008783,
  0.9586550990404803,
  0.9633846254394740,
  0.9678371941389582,
  0.9720101439151101,
  0.9759009806706322,
  0.9795073789268500,
  0.9828271832159826,
  0.9858584093735683,
  0.9885992457319537,
  0.9910480542178592,
  0.9932033713622931,
  0.9950639092458672,
  0.9966285564501065,
  0.9978963792674906,
  0.9988666243127571,
  0.9995387299886880,
  0.9999124439735659
};

double GL_weight[] = {
  0.0002246904801461,
  0.0005229063396701,
  0.0008212515093345,
  0.0011191442154813,
  0.0014163757357290,
  0.0017127630204551,
  0.0020081274918693,
  0.0023022921283514,
  0.0025950809163382,
  0.0028863187714329,
  0.0031758315808536,
  0.0034634462834494,
  0.0037489909628174,
  0.0040322949452430,
  0.0043131888993084,
  0.0045915049358304,
  0.0048670767075034,
  0.0051397395079161,
  0.0054093303697515,
  0.0056756881620402,
  0.0059386536863701,
  0.0061980697719754,
  0.0064537813696337,
  0.0067056356443082,
  0.0069534820664760,
  0.0071971725020834,
  0.0074365613010736,
  0.0076715053844326,
  0.0079018643296997,
  0.0081275004548926,
  0.0083482789007946,
  0.0085640677115557,
  0.0087747379135588,
  0.0089801635925043,
  0.0091802219686657,
  0.0093747934702723,
  0.0095637618049755,
  0.0097470140293533,
  0.0099244406164154,
  0.0100959355210650,
  0.0102613962434801,
  0.0104207238903756,
  0.0105738232341107,
  0.0107206027696042,
  0.0108609747690261,
  0.0109948553342302,
  0.0111221644468999,
  0.0112428260163725,
  0.0113567679251182,
  0.0114639220718434,
  0.0115642244121935,
  0.0116576149970314,
  0.0117440380082680,
  0.0118234417922238,
  0.0118957788905017,
  0.0119610060683517,
  0.0120190843405120,
  0.0120699789945096,
  0.0121136596114076,
  0.0121501000839859,
  0.0121792786323453,
  0.0122011778169248,
  0.0122157845489250,
  0.0122230900981313,
  0.0122230900981313,
  0.0122157845489250,
  0.0122011778169248,
  0.0121792786323453,
  0.0121501000839859,
  0.0121136596114076,
  0.0120699789945096,
  0.0120190843405120,
  0.0119610060683517,
  0.0118957788905017,
  0.0118234417922238,
  0.0117440380082680,
  0.0116576149970314,
  0.0115642244121935,
  0.0114639220718434,
  0.0113567679251182,
  0.0112428260163725,
  0.0111221644468999,
  0.0109948553342302,
  0.0108609747690261,
  0.0107206027696042,
  0.0105738232341107,
  0.0104207238903756,
  0.0102613962434801,
  0.0100959355210650,
  0.0099244406164154,
  0.0097470140293533,
  0.0095637618049755,
  0.0093747934702723,
  0.0091802219686657,
  0.0089801635925043,
  0.0087747379135588,
  0.0085640677115557,
  0.0083482789007946,
  0.0081275004548926,
  0.0079018643296997,
  0.0076715053844326,
  0.0074365613010736,
  0.0071971725020834,
  0.0069534820664760,
  0.0067056356443082,
  0.0064537813696337,
  0.0061980697719754,
  0.0059386536863701,
  0.0056756881620402,
  0.0054093303697515,
  0.0051397395079161,
  0.0048670767075034,
  0.0045915049358304,
  0.0043131888993084,
  0.0040322949452430,
  0.0037489909628174,
  0.0034634462834494,
  0.0031758315808536,
  0.0028863187714329,
  0.0025950809163382,
  0.0023022921283514,
  0.0020081274918693,
  0.0017127630204551,
  0.0014163757357290,
  0.0011191442154813,
  0.0008212515093345,
  0.0005229063396701,
  0.0002246904801461
};

double rk2_leg ( double t1, double t2, double x, int n )
/*
  RK2_LEG advances the value of X(T) using a Runge-Kutta method.
*/
{
  double f,h,k1,k2,snn1,t;
  int j,m=10;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );

  t = t1;

  for ( j = 0; j < m; j++ )
  {
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;

    t = t + h;

    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );   
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
}

double ts_mult ( double *u, double h, int n )
/*
    TS_MULT evaluates a polynomials
*/
{
  double hk,ts;
  int k;
  
  ts = 0.0;
  hk = 1.0;
  for ( k = 1; k<= n; k++ )
  {
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}


void legendre_compute_glr0 ( int n, double *p, double *pp )

/*
    LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
*/
{
  double dk;
  int k;
  double pm1,pm2,ppm1,ppm2;

  pm2 = 0.0;
  pm1 = 1.0;
  ppm2 = 0.0;
  ppm1 = 0.0;

  for ( k = 0; k < n; k++ )
  {
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}

void legendre_compute_glr1 ( int n, double *x, double *ders )
/*
    LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
*/
{
  double dk,dn,h;
  int j,k,l,m=30,n2,s;

  double *u;
  double *up;
  double xp;

  if ( n % 2 == 1 )
  {
    n2 = ( n - 1 ) / 2;
    s = 1;
  }
  else
  {
    n2 = n / 2;
    s = 0;
  }

  u  = ( double * ) _mm_malloc ( ( m + 2 ) * sizeof ( double ) ,32);
  up = ( double * ) _mm_malloc ( ( m + 1 ) * sizeof ( double ) ,32);

  dn = ( double ) n;

  for ( j = n2; j < n - 1; j++ )
  {
    xp = x[j];

    h = rk2_leg ( pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = ders[j];

    up[0] = 0.0;
    up[1] = u[2];

    for ( k = 0; k <= m - 2; k++ )
    {
      dk = ( double ) k;

      u[k+3] = 
      ( 
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 )
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }

    for ( l = 0; l < 5; l++ )
    { 
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    ders[j+1] = ts_mult ( up, h, m-1 );
  }

  free ( u );
  free ( up );

  for ( k = 0; k < n2 + s; k++ )
  {
    x[k] = - x[n-k-1];
    ders[k] = ders[n-k-1];
  }
  return;
}

void legendre_compute_glr2 ( double pn0, int n, double *x1,  double *d1 )
/*
    LEGENDRE_COMPUTE_GLR2 finds the first real root.
*/
{
  double dk,dn,t;
  int k,l,m=30;

  double *u;
  double *up;

  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u  = ( double * ) _mm_malloc ( ( m + 2 ) * sizeof ( double ) ,32);
  up = ( double * ) _mm_malloc ( ( m + 1 ) * sizeof ( double ) ,32);

  dn = ( double ) n;
/*
  U[0] and UP[0] are never used.
  U[M+1] is set, but not used, and UP[M] is set and not used.
  What gives?
*/
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;
 
  for ( k = 0; k <= m - 2; k = k + 2 )
  {
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / ( dk + 1.0 ) / ( dk + 2.0 );
 
    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }
  
  for ( l = 0; l < 5; l++ )
  {
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  free ( u );
  free ( up) ;

  return;
}


void legendre_compute_glr ( int n, double x[], double w[] )
/*
    LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
*/
{
  int i;
  double p=0;
  double pp=0;
  double w_sum;
/*
  Get the value and derivative of the N-th Legendre polynomial at 0.
*/
  legendre_compute_glr0 ( n, &p, &pp );
/*
  Either zero is a root, or we have to call a function to find the first root.
*/  
  if ( n % 2 == 1 )
  {
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }
  else
  {
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
/*
  Get the complete set of roots and derivatives.
*/
  legendre_compute_glr1 ( n, x, w );
/*
  Compute the weights.
*/
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}

void rescale ( double a, double b, int n, double x[], double w[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}

void glwt_fast( int n, double a, double b , double *x, double* w)
/*
    LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
*/ 
{
 /*
  Compute the rule.
*/
  legendre_compute_glr ( n, x, w );
/*
  Rescale the rule to [A,B].
*/
  rescale ( a, b, n, x, w );
}

void glwt(int n, double a, double b, double *x, double *w)
{
  int i,k,N,N1,N2;
  double h,d,H;
  N  = n - 1;
  N1 = N + 1; 
  N2 = N + 2;
  
  double *xu = (double*) _mm_malloc(N1*sizeof(double),32);
  double *y  = (double*) _mm_malloc(N1*sizeof(double),32);
  double *L  = (double*) _mm_malloc(N1*sizeof(double)*3,32);
  double *Lp = (double*) _mm_malloc(N1*sizeof(double),32);
  // Initial guess
  h = 2./(double) (N1-1.);
  for (i=0;i<N1;i++)
    {
      xu[i]   = -1. + i*h;
      y[i] = cos( (2.*i+1.)*pi/(2.*N+2.) )+ 
	(0.27)/(double) N1*sin(pi*xu[i]*N/(double) N2);
    }
  d = 2;
  while (d>EPSILON)
    {
      for(i=0;i<N1;i++)
	{
	  L[i] = 1.;
	  L[N1+i] = y[i];
	}
      for (k=1;k<N1-1;k++)
	for(i=0;i<N1;i++)
	  {
	    L[2*N1+i] = ( (2.*k+1.)*y[i]*L[N1+i]-k*L[i] )/(double) (k+1.);
	    L[i]    = L[N1+i];
	    L[N1+i] = L[2*N1+i];
	  }
      // last iteration remained
      k = N1-1;
      for(i=0;i<N1;i++)
	{
	  L[2*N1+i] = ( (2.*k+1.)*y[i]*L[N1+i]-k*L[i] )/(double) (k+1.);
	}
      d = 0;
      for(i=0;i<N1;i++)
	{
	  Lp[i] = N2*( L[N1+i]-y[i]*L[2*N1+i]  )/( 1-y[i]*y[i]  );
	  H = L[2*N1+i]/Lp[i];
	  d = MAX(d,fabs(H));
	  y[i] -= H;
	}
    }
 
  // Linear map from [-1,1] to [a,b]
  double c = N2*N2/(double) (N1*N1);
  for (i=0;i<N1;i++)
    {
      x[N1-i-1] = ( a*(1.-y[i])+b*(1.+y[i]) )/2.;
      w[N1-i-1] = (b-a)/( (1.-y[i]*y[i])*Lp[i]*Lp[i] )*c;
    }
}

void trapz_wt(int n, double *x, double *w)
{
  int i;
  double h = 1./(double) (n-1);
 
  for (i=0;i<n;i++)
    {
      x[i] = i*h;
      w[i] = h;
    }
  w[0  ] /=2.;
  w[n-1] /=2.;
}

static inline double K0(double a, double b, double x)
{
  if(x==0)
    return 0.;
  else
    return 1./x*exp(-b*x-a/x);
}

double IncompBesselK0_int(double a, double b, int n,double *x, double*w)
{
#ifdef __AVX__
  int i;
  double s[4] __attribute__((aligned(32)));
  
  __m256d S,X,C,A,B,AX,BX,ONE,InvX,ExpC,W;
  S = _mm256_setzero_pd();

  for (i=0;i<n;i+=4)
    {
      X    = _mm256_load_pd( x + i );
      B    = _mm256_set1_pd( -b );
      A    = _mm256_set1_pd( -a );
      AX   = _mm256_div_pd( A, X  );
      BX   = _mm256_mul_pd( B, X  );
      C    = _mm256_add_pd( AX, BX  );
      ONE  = _mm256_set1_pd( 1.0 );
      InvX = _mm256_div_pd(ONE , X);
      ExpC = se_mm256_exp_pd( C );
      W    = _mm256_load_pd( w + i );
      S    = _mm256_add_pd(S, _mm256_mul_pd( W, _mm256_mul_pd(InvX, ExpC) ));
    }
  _mm256_store_pd(s, S);
  return s[0]+s[1]+s[2]+s[3];
#elif defined __SSE4_2__
  int i;
  double s[2] __attribute__((aligned(16)));
  
  __m128d S,X,C,A,B,AX,BX,ONE,InvX,ExpC,W;
  S = _mm_setzero_pd();

  for (i=0;i<n;i+=2)
    {
      X    = _mm_load_pd( x + i );
      B    = _mm_set1_pd( -b );
      A    = _mm_set1_pd( -a );
      AX   = _mm_div_pd( A, X  );
      BX   = _mm_mul_pd( B, X  );
      C    = _mm_add_pd( AX, BX  );
      ONE  = _mm_set1_pd( 1.0 );
      InvX = _mm_div_pd(ONE , X);
      ExpC = se_mm_exp_pd( C );
      W    = _mm_load_pd( w + i );
      S    = _mm_add_pd(S, _mm_mul_pd( W, _mm_mul_pd(InvX, ExpC  ) ));
    }
  _mm_store_pd(s, S);
  return s[0]+s[1];
#else
  int i;
  double s=0.,f;
  for (i=0;i<n;i++)
    {
      f = K0(a,b,x[i]);
      s+= w[i]*f;
    }
  return s;
#endif
}


double IncompBesselK0_int_inf(double a, double b, int n)
{
  int i;
  double eeps = 1.e-15;
  double aa = log(1./eeps)/a;
  double bb = b/a;
  
  unsigned long int tmax = aa/2+sqrt(aa*aa/4-bb);
  double h = (double) (tmax-1.)/n;
  
  double *tvec = (double*) _mm_malloc(sizeof(double)*(n+1),32);
  double *wv   = (double*) _mm_malloc(sizeof(double)*(n+1),32);
  for (i=0; i<=n;i++)
    {
      tvec[i] = 1.+i*h;
      wv[i]   = 1.;
    }
  wv[0] = .5; wv[n] = .5;
  
  double  corrv[5] ={-245./1440., 
		     462./1440.,
		     -336./1440.,
		     146./1440.,
		     -27./1440.};
  
  for (i=0;i<5;i++)
    {
      wv[i] += corrv[i];
      wv[n-4+i]+=corrv[5-i];
    }
  
  double K0val=0.;
  for (i=0;i<=n;i++)
    K0val += h*wv[i]*exp(-a*tvec[i]-b/tvec[i])/tvec[i];

  return K0val;

}


/* Incomplete Bessel with Adaptive Simpson rule*/
double f(double x, double bes_a, double bes_b, int der)
{
  if(x==0)
    return 0;
  if(der==1)
     return exp(-bes_a/x-bes_b*x);
  else if(der==0)
     return exp(-bes_a/x-bes_b*x)/x;
  return 0;
}

double
simpson(double a, double b, double bes_a, double bes_b, int der)
{
  double c = (a+b)/2.;
  return (b-a)/6.*(f(a, bes_a, bes_b, der)+4.*f(c, bes_a, bes_b, der)+f(b, bes_a, bes_b, der));
}

double
errest(double a, double b, double *val, double bes_a, double bes_b, int der)
{
  double c = (a+b)/2.;
  double I = simpson(a,b, bes_a, bes_b, der);
  *val      = simpson(a,c, bes_a, bes_b, der) + simpson(c,b, bes_a, bes_b, der);
  return fabs(I- (*val) );
}

double
IncompBesselK0_Simpson(double tol, int *cnt, double bes_a, double bes_b, int der)
{

  int HEAD, n;
  double res = 0.0, val;
  double a, b;
  /* FIXME: size of stack should be dynamic. One possibility is 
   * is to use vectors. Otherwise one can keep only the last value 
   * of the interval and if OK reset it to c2 again! (expensive)
   */
  double *ST = (double*) malloc(STACK_SIZE*sizeof(double));

  // First value in the stack is upper bound of the integral
  ST[0] = 1.;
  
  if( 15.*errest(0.,1.,&val, bes_a, bes_b,der)<tol)
    return val;
  a = 0.;
  b = 1.;
  
  HEAD  = 0;
  
  n=0;  // to count the number of splits
  for(;;)
    {
      if(15.*errest(a,b,&val, bes_a, bes_b,der)>tol)
        {
          if ( HEAD>=(STACK_SIZE-1) )
            {
              printf("Increase the stack size!\n");
              exit(EXIT_FAILURE);
            }
          ST[++HEAD] = (a+b)/2.;
          b          = ST[HEAD];
          n++;
        }
      else
        {
          res += val;
          // if a==b or HEAD is pointing the the first element we break
          if(HEAD == 0)
            break;
          a    = ST[HEAD];
          b    = ST[--HEAD];
        }
    }
  free(ST);
  *cnt = n;
  return res;
}

/* kernels for computing the modified bessel function of the 
 * second kind and its derivative */
double bessel_f(double x, void * p)
{
  gsl_params * params = (gsl_params *)p;
  double a = (params->a);
  double b = (params->b);
  return exp(-a/x-b*x)/x;
}

double bessel_f_der(double x, void * p)
{
  gsl_params * params = (gsl_params *)p;
  double a = (params->a);
  double b = (params->b);
  return exp(-a/x-b*x);
}

double call_gsl_bessel_integrator(double a, double b, 
				  gsl_integration_workspace *w,
				  int der)
{
  double result, error;
  gsl_function F;
  if(der==0)
    F.function = &bessel_f;
  else if(der==1)
    F.function = &bessel_f_der;

  if(a==0){
    result = 1.e308;
    return result;
  }

  if(a<b && 0){
    gsl_params  params = {b,a};
    F.params = &params;
    gsl_integration_qags (&F, 0, 1, 0, 1e-12, 10000,
     			  w, &result, &error);
    double z = 2.*gsl_sf_bessel_K0(2.*sqrt(a*b));
    result = z-result;
  }
  else{
    gsl_params  params = {a,b};
    F.params = &params;
    gsl_integration_qags (&F, 0, 1, 0, 1e-12, 10000, 
			  w, &result, &error); 
  }
  
  return result;   
}


  
double func(double x, double a, double b)
{
  if(x==0)
    return 0;
  else
    return exp(-a/x-b*x)/x;
}

double computeK0(double a, double b) {
    // integration over [0,1] with "panels" panels
  double s=0; 
  int panels = 4;
  double x0;
  double h = 1./panels;
  for (int p=0;p<panels;p++) {
    x0 = 0+p*h;
#ifdef _OPENMP
#pragma omp parallel for  reduction(+:s)
#endif
    for (int i=0; i<128;i++) {
      double x = x0 + h*GL_value[i];
      s += GL_weight[i]*func(x,a,b)*h;
    }
  }
  return s;
}
