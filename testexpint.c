#include <stdio.h>
#include "math.h"
#include "immintrin.h"
#include "math_x86.h"
#include "expint.h"



__m256d _mm256_ein_pd(__m256d Z)
{

  


  return Z;
}










int main()
{
  int i;
  double x[] = {.7, 1.3, 2.7, 6.};

  for (i=0; i<4; i++)
    printf("%f ",expint_log_euler(x[i]));
  printf("\n");

  double s[4] __attribute__((aligned(32)));
  __m256d Z;
  Z = _mm256_set_pd(x[3],x[2],x[1],x[0]);
  Z = _mm256_ein_pd(Z);
  _mm256_store_pd(s,Z);
  for (i=0; i<4; i++)
    printf("%f ",expint_log_euler(s[i]));
  printf("\n");
  return 0;
}
