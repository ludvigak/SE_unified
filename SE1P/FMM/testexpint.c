#include <stdio.h>
#include "math.h"
#include "immintrin.h"
#include "math_x86.h"
#include "expint.h"
#include "stdbool.h"
#include "mex_FMM.h"

int main()
{
  int i;
  unsigned int N = 3000;
  double *x = (double*) _mm_malloc(N*sizeof(double),32);
  double *ymat = (double*) _mm_malloc(N*sizeof(double),32);
  double *yavx = (double*) _mm_malloc(N*sizeof(double),32);

  __m256d Z;
  srand(time(NULL));
  for (i=0; i<N; i++)
    x[i] = 10.0*rand()/RAND_MAX;

  if(N<20)
    for (i=0; i<N; i++)
      {
	printf("%f %f\n",expint_log_euler(x[i]),expints(1,x[i])+log(x[i])+EULER);
      }


  /*matlab ein(x)*/
  double start = gettime();
  for (i=0; i<N; i++)
    ymat[i] = expint_log_euler(x[i]);
  printf("MATLAB %f\n",gettime()-start);

  start = gettime();
  /*AVX ein(x)*/
  __m128d Z0;
  for (i=0; i<N; i+=8)
    {
      int i0=0;
      Z = _mm256_set_pd(x[i+3+i0],x[i+2+i0],x[i+1+i0],x[i+i0]);
      Z = _mm256_ein_pd(Z);
      _mm256_store_pd(yavx+i+i0,Z);
      i0+=4;
      Z = _mm256_set_pd(x[i+3+i0],x[i+2+i0],x[i+1+i0],x[i+i0]);
      Z = _mm256_ein_pd(Z);
      _mm256_store_pd(yavx+i+i0,Z);
    }
  /* remainder loop*/
  for (int j=i; j<N; j++)
    yavx[j] = expint_log_euler(x[j]);
  printf("AVX    %f\n",gettime()-start);

  double s0 = 0, s1=0;
  for (i=0; i<N; i++)
    {
      s0+=ymat[i]*ymat[i];
      s1+= (ymat[i]-yavx[i])*(ymat[i]-yavx[i]);
    }
  printf("%g\n",sqrt(s1/s0));
  return 0;
}
