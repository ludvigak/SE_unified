#include "x86intrin.h"
#include "string.h"
#include "malloc.h"
#include "SE_fgg.h"
#include "SE_fkg.h"


// Kaiser window function gridding
// Davoud Saffar Shamshirgar, davoudss@kth.se


// =============================================================================
// Core SE FGG routines ========================================================

void SE_FKG_allocate_workspace(SE_FGG_work* work, const SE_FGG_params* params, 
			       int allocate_fgg_expa)
{
    int numel = SE_prod3(params->npdims);
    work->H = SE_FGG_MALLOC(numel*sizeof(double));
  
    SE_fp_set_zero(work->H, numel);
    
    work->free_zs = false;
    work->zs = NULL;
  
    if(allocate_fgg_expa)
    {
	numel = (params->N)*(params->P);
	work->zx = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->zy = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->zz = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->idx= (int*) SE_FGG_MALLOC(params->N*sizeof(int));
    }
    else
    {
	work->zx=NULL;
	work->zy=NULL;
	work->zz=NULL;
	work->idx=NULL;
    }
    work->free_fgg_expa=allocate_fgg_expa;
}

// -----------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ----------------KAISER KERNEL ---------------------------------------------
// ---------------------------------------------------------------------------

static inline double
kaiser(double x, double ow2, double beta) {
  double t = sqrt(1. - x*x*ow2);
  return exp(beta*(t-1));
}

static
int kaiser_expansion_3p(const double x[3], const double q,
			const SE_FGG_params* params,
			double z2_0[P_MAX],
			double z2_1[P_MAX],
			double z2_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double one_h = 1/h;
    const double w = params->P/2.;
    const double ow2  =1./(w*w);
    const double beta = params->beta;
    double t0[3];

    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p)) {
      for(int j=0; j<3; j++)
	{
	  idx = (int) round(x[j]*one_h);
	  idx_from[j] = idx - p_half;
	  t0[j] = x[j]*one_h-idx_from[j];
	}
    }
    else {
      for(int j=0; j<3; j++)
	{
	  idx = (int) floor(x[j]*one_h);
	  idx_from[j] = idx - (p_half-1);
	  t0[j] = x[j]*one_h-idx_from[j];
	}
    }

    // compute second factor by induction
    for(int i=0; i<p; i++) {
      z2_0[i] = kaiser(t0[0]-i,ow2,beta);
      z2_1[i] = kaiser(t0[1]-i,ow2,beta);
      z2_2[i] = kaiser(t0[2]-i,ow2,beta);
    }

    // save some flops by multiplying one vector with q
    for(int i=0; i<p; i++)
      z2_0[i] *= q;
    
    return __IDX3_RMAJ(idx_from[0]+p_half,
                       idx_from[1]+p_half,
                       idx_from[2]+p_half,
                       params->npdims[1], params->npdims[2]);
}

// -----------------------------------------------------------------------------
void SE_FKG_expand_all(SE_FGG_work* work, 
			      const SE_state* st, 
			      const SE_FGG_params* params)
{
    double xn[3] MEM_ALIGNED;
    const int N = params->N;
    const int P = params->P;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
      // compute index and expansion vectors
      xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
      
      *(work->idx+n) = kaiser_expansion_3p(xn,1,params, 
					   work->zx+n*P, 
					   work->zy+n*P, 
					   work->zz+n*P);
    }
}


// -----------------------------------------------------------------------------
void SE_FKG_int(double* restrict phi,  
		const SE_FGG_work* work, 
		const SE_state* st, 
		const SE_FGG_params* params)
{
    double z2_0[P_MAX] MEM_ALIGNED;
    double z2_1[P_MAX] MEM_ALIGNED;
    double z2_2[P_MAX] MEM_ALIGNED;

    // unpack params
    const double* restrict H = work->H;
    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    double xm[3];
    int i,j,k,idx;
    double phi_m, cij;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++) {
      xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N];
      
      idx = kaiser_expansion_3p(xm, 1, params, z2_0, z2_1, z2_2);
      
      phi_m = 0;
      
      for(i = 0; i<p; i++) {
	for(j = 0; j<p; j++) {
	  cij = z2_0[i]*z2_1[j];
	  for(k = 0; k<p; k++) {
	    phi_m += H[idx]*z2_2[k]*cij;
	    idx++;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
      phi[m] = (h*h*h)*phi_m;
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_int_split(double* restrict phi,  
		      const SE_FGG_work* work, 
		      const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,k,idx,idx_zz;
    double phi_m, cij;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	phi_m = 0;

	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij = zx[m*p+i]*zy[m*p+j];
		idx_zz=m*p;
		for(k = 0; k<p; k++)
		{
		    phi_m += H[idx]*zz[idx_zz]*cij;
		    idx++; idx_zz++;
		}
		idx += incrj;
	    }
	    idx += incri;
	}
	phi[m] = h3*phi_m;
    }
}

void SE_FKG_int_split_SSE(double* restrict phi,  
			  const SE_FGG_work* work, 
			  const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,k,idx,idx_zz;
    double s[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rC, rP;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];	
	rP=_mm_setzero_pd();
	int mp = m*p;
	for(i = 0; i<p; i++)
	  {
	    for(j = 0; j<p; j++)
	      {
		rC = _mm_set1_pd( zx[mp+i]*zy[mp+j]);
		idx_zz=mp;
		for(k = 0; k<p; k+=2)
		  {
		    rH0  = _mm_loadu_pd( H+idx );
		    rZZ0 = _mm_load_pd( zz + idx_zz);
#ifdef AVX_FMA
		    rP = _mm_fmadd_pd(rH0,_mm_mul_pd(rZZ0,rC),rP);
#else
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rC)));
#endif

		    idx+=2; 
		    idx_zz+=2;
		  }
		idx += incrj;
	      }
	    idx += incri;
	  }
	_mm_store_pd(s,rP);
	phi[m] = h3*(s[0]+s[1]);
    }
}

void SE_FKG_int_split_SSE_P6(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx,idx_zz;
    double s[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rC, rP, rT;
    __m128d rH1, rZZ1;
    __m128d rH2, rZZ2;

    const int incrj = params->npdims[2]-6;
    const int incri = params->npdims[2]*(params->npdims[1]-6);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];	
	rP=_mm_setzero_pd();
	int sixm = 6*m;
	for(i = 0; i<6; i++)
	  {
	    for(j = 0; j<6; j++)
	      {
		rC = _mm_set1_pd( zx[sixm+i]*zy[sixm+j]);
		idx_zz= sixm;

		rH0  = _mm_loadu_pd( H+idx    );
		rH1  = _mm_loadu_pd( H+idx + 2);
		rH2  = _mm_loadu_pd( H+idx + 4);
		
		rZZ0 = _mm_load_pd( zz + idx_zz    );
		rZZ1 = _mm_load_pd( zz + idx_zz + 2);
		rZZ2 = _mm_load_pd( zz + idx_zz + 4);

#ifdef AVX_FMA
		rP = _mm_fmadd_pd(rH0,_mm_mul_pd(rZZ0,rC),rP);
		rP = _mm_fmadd_pd(rH1,_mm_mul_pd(rZZ1,rC),rP);
		rP = _mm_fmadd_pd(rH2,_mm_mul_pd(rZZ2,rC),rP);
#else
		rT = _mm_add_pd(_mm_mul_pd(rH0,rZZ0),_mm_mul_pd(rH1,rZZ1));
		rT = _mm_mul_pd(_mm_add_pd(rT,_mm_mul_pd(rH2,rZZ2)),rC);

		rP = _mm_add_pd(rP,rT);
#endif		
		
		idx += incrj + 6;
	      }
	    idx += incri;
	  }
	_mm_store_pd(s,rP);
	phi[m] = h3*(s[0]+s[1]);
    }
}

void SE_FKG_int_split_SSE_P8(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=8 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx;
    double s[2] MEM_ALIGNED;

    // hold entire zz vector
    __m128d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m128d rC, rP;
    __m128d rH0, rH1, rH2, rH3; 

    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	rP=_mm_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*8     );
	rZZ1 = _mm_load_pd(zz + m*8 + 2 );
	rZZ2 = _mm_load_pd(zz + m*8 + 4 );
	rZZ3 = _mm_load_pd(zz + m*8 + 6 );

	for(i = 0; i<8; i++)
	  {
	    for(j = 0; j<8; j++)
	      {
		rC = _mm_set1_pd( zx[m*8+i]*zy[m*8+j]);

		rH0  = _mm_loadu_pd( H+idx    );
		rH1  = _mm_loadu_pd( H+idx + 2);
		rH2  = _mm_loadu_pd( H+idx + 4);
		rH3  = _mm_loadu_pd( H+idx + 6);
#ifdef AVX_FMA
		rP = _mm_fmadd_pd(rH0,_mm_mul_pd(rZZ0,rC),rP);
		rP = _mm_fmadd_pd(rH1,_mm_mul_pd(rZZ1,rC),rP);
		rP = _mm_fmadd_pd(rH2,_mm_mul_pd(rZZ2,rC),rP);
		rP = _mm_fmadd_pd(rH3,_mm_mul_pd(rZZ3,rC),rP);
#else		
		rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rC)));
#endif
		
		idx += incrj + 8;
	      }
		idx += incri;
	  }
	_mm_store_pd(s,rP);
	phi[m] = h3*(s[0]+s[1]);
    }
}

void SE_FKG_int_split_SSE_P16(double* restrict phi,  
			      const SE_FGG_work* work, 
			      const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx;
    double s[2] MEM_ALIGNED;

    // hold entire zz vector
    __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
    __m128d rC, rP;
    __m128d rH0, rH1, rH2, rH3;

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

	rP=_mm_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*16     );
	rZZ1 = _mm_load_pd(zz + m*16 + 2 );
	rZZ2 = _mm_load_pd(zz + m*16 + 4 );
	rZZ3 = _mm_load_pd(zz + m*16 + 6 );
	rZZ4 = _mm_load_pd(zz + m*16 + 8 );
	rZZ5 = _mm_load_pd(zz + m*16 + 10);
	rZZ6 = _mm_load_pd(zz + m*16 + 12);
	rZZ7 = _mm_load_pd(zz + m*16 + 14);

	for(i = 0; i<16; i++)
	  {
	    for(j = 0; j<16; j++)
	      {
		rC = _mm_set1_pd( zx[m*16+i]*zy[m*16+j]);

		rH0  = _mm_loadu_pd( H+idx     );
		rH1  = _mm_loadu_pd( H+idx + 2 );
		rH2  = _mm_loadu_pd( H+idx + 4 );		
		rH3  = _mm_loadu_pd( H+idx + 6 );
#ifdef AVX_FMA	
		rP = _mm_fmadd_pd(rH0,_mm_mul_pd(rZZ0,rC),rP);
		rP = _mm_fmadd_pd(rH1,_mm_mul_pd(rZZ1,rC),rP);
		rP = _mm_fmadd_pd(rH2,_mm_mul_pd(rZZ2,rC),rP);
		rP = _mm_fmadd_pd(rH3,_mm_mul_pd(rZZ3,rC),rP);
#else		
		rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rC)));
#endif
		rH0  = _mm_loadu_pd( H+idx + 8 );
		rH1  = _mm_loadu_pd( H+idx + 10);
		rH2  = _mm_loadu_pd( H+idx + 12);
		rH3  = _mm_loadu_pd( H+idx + 14);
#ifdef AVX_FMA	
		rP = _mm_fmadd_pd(rH0,_mm_mul_pd(rZZ4,rC),rP);
		rP = _mm_fmadd_pd(rH1,_mm_mul_pd(rZZ5,rC),rP);
		rP = _mm_fmadd_pd(rH2,_mm_mul_pd(rZZ6,rC),rP);
		rP = _mm_fmadd_pd(rH3,_mm_mul_pd(rZZ7,rC),rP);
#else		
		rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(rZZ4,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(rZZ5,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(rZZ6,rC)));
		rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(rZZ7,rC)));
#endif
		
		idx += incrj + 16;
	      }
	    idx += incri;
	  }
	_mm_store_pd(s,rP);
	phi[m] = h3*(s[0]+s[1]);
    }
}

void SE_FKG_int_split_AVX_P4(double* restrict phi,
                             const SE_FGG_work* work,
                             const SE_FGG_params* params)
{  
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    // ASSUME P=8 const int p = params->P; //
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx;
    double s[4] MEM_ALIGNED;
    
    // hold entire zz vector
    __m256d rZZ0;
    __m256d rC, rP;
    __m256d rH0;
    
    const int incrj = params->npdims[2]-4;
    const int incri = params->npdims[2]*(params->npdims[1]-4);
    
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        rP = _mm256_setzero_pd();
	int m4 = m*4;
        // hoist load of ZZ vector //
        rZZ0 = _mm256_load_pd(zz + m4     );

	for(i = 0; i<4; i++)
	  {
	    for(j = 0; j<4; j++)
	      {
		rC = _mm256_set1_pd( zx[m4+i]*zy[m4+j]);

		rH0  = _mm256_loadu_pd( H+idx    );

#ifdef AVX_FMA
		rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(rZZ0,rC),rP);
#else
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(rZZ0,rC)));
#endif
		idx += incrj + 4;
	      }
	    idx += incri;

        }
        _mm256_store_pd(s,rP);

	phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
    }
}

void SE_FKG_int_split_AVX_P6(double* restrict phi,
                             const SE_FGG_work* work,
                             const SE_FGG_params* params)
{  
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    // ASSUME P=8 const int p = params->P; //
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx;
    double s[4] MEM_ALIGNED;
    
    // hold entire zz vector
    __m256d rZZ0, rC0, rP0, rH0;
    __m128d rZZ1, rH1, rC1, rP1; 
    
    const int incrj = params->npdims[2]-6;
    const int incri = params->npdims[2]*(params->npdims[1]-6);
    
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        rP0 = _mm256_setzero_pd();
	rP1 = _mm_setzero_pd();

        // hoist load of ZZ vector //
        rZZ0 = _mm256_loadu_pd(zz + m*6     );
        rZZ1 = _mm_loadu_pd(zz + m*6 + 4 );
	int sixm = 6*m;
	for(i = 0; i<6; i++)
	  {
	    double tmp = zx[sixm+i];
	    for(j = 0; j<6; j++)
	      {
		rC0 = _mm256_set1_pd( tmp*zy[sixm+j]);
		rC1 = _mm_set1_pd( tmp*zy[sixm+j]);

		rH0  = _mm256_loadu_pd( H+idx    );
		rH1  = _mm_loadu_pd(    H+idx + 4);
#ifdef AVX_FMA
		rP0 = _mm256_fmadd_pd(rH0,_mm256_mul_pd(rZZ0,rC0),rP0);
		rP1 = _mm_fmadd_pd(   rH1,_mm_mul_pd(   rZZ1,rC1),rP1);
#else
		rP0 = _mm256_add_pd(rP0,_mm256_mul_pd(rH0,_mm256_mul_pd(rZZ0,rC0)));
		rP1 = _mm_add_pd(rP1,_mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rC1)));
#endif

		idx += incrj + 6;
	      }
	    idx += incri;

        }
        _mm256_store_pd(s,rP0);
	phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
	
        _mm_store_pd(s,rP1);
	phi[m] += h3*(s[0]+s[1]);

    }
}

void SE_FKG_int_split_AVX_P8(double* restrict phi,
                             const SE_FGG_work* work,
                             const SE_FGG_params* params)
{  
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    // ASSUME P=8 const int p = params->P; //
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx;
    double s[4] MEM_ALIGNED;
    
    // hold entire zz vector
    __m256d rZZ0, rZZ1;
    __m256d rC, rP;
    __m256d rH0, rH1;
    
    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);
    
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        rP = _mm256_setzero_pd();

        // hoist load of ZZ vector //
        rZZ0 = _mm256_load_pd(zz + m*8     );
        rZZ1 = _mm256_load_pd(zz + m*8 + 4 );

	for(i = 0; i<8; i++)
	  {
	    for(j = 0; j<8; j++)
	      {
		rC = _mm256_set1_pd( zx[m*8+i]*zy[m*8+j]);

		rH0  = _mm256_loadu_pd( H+idx    );
		rH1  = _mm256_loadu_pd( H+idx + 4);

#ifdef AVX_FMA
		rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(rZZ0,rC),rP);
		rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(rZZ1,rC),rP);
#else
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(rZZ0,rC)));
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(rZZ1,rC)));
#endif
		idx += incrj + 8;
	      }
	    idx += incri;

        }
        _mm256_store_pd(s,rP);

	phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_int_split_AVX_P16(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;    

    int i,j,idx;
    double s[4] MEM_ALIGNED;

    // hold entire zz vector
    __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m256d rC, rP;
    __m256d rH0, rH1, rH2, rH3; 

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	rP=_mm256_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm256_load_pd(zz + m*16     );
	rZZ1 = _mm256_load_pd(zz + m*16 + 4 );
	rZZ2 = _mm256_load_pd(zz + m*16 + 8 );
	rZZ3 = _mm256_load_pd(zz + m*16 + 12);

	for(i = 0; i<16; i++)
	  {
	    for(j = 0; j<16; j++)
	      {
		rC = _mm256_set1_pd( zx[m*16+i]*zy[m*16+j]);

		rH0  = _mm256_loadu_pd( H+idx     );
		rH1  = _mm256_loadu_pd( H+idx + 4 );
		rH2  = _mm256_loadu_pd( H+idx + 8 );
		rH3  = _mm256_loadu_pd( H+idx + 12);

#ifdef AVX_FMA
		rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(rZZ0,rC),rP);
		rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(rZZ1,rC),rP);
		rP = _mm256_fmadd_pd(rH2,_mm256_mul_pd(rZZ2,rC),rP);
		rP = _mm256_fmadd_pd(rH3,_mm256_mul_pd(rZZ3,rC),rP);
#else		    
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(rZZ0,rC)));
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(rZZ1,rC)));
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(rZZ2,rC)));
		rP = _mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(rZZ3,rC)));
#endif
		idx += incrj + 16;
	      }
	    idx += incri;
	  }
	_mm256_store_pd(s,rP);
	phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
    }
}


void SE_FKG_int_split_SSE_dispatch(double* restrict phi,  
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    __DISPATCHER_MSG("[FKG INT SSE] SSE Disabled\n");
    SE_FKG_int_split(phi, work, params);
    return;
#endif

    // if P is odd, or if either increment is odd, fall back on vanilla
    if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
	__DISPATCHER_MSG("[FKG INT SSE] SSE Abort (PARAMS)\n");
	SE_FKG_int_split(phi, work, params);
	return;
    }

#if 0
    if(	( (unsigned long) work->zx)%16 != 0 || 
	( (unsigned long) work->zy)%16 != 0 ||
	( (unsigned long) work->zz)%16 != 0 )
    {
	__DISPATCHER_MSG("[FKG INT SSE] SSE Abort (DATA)\n");
	SE_FKG_int_split(phi, work, params);
	return;
    }
#endif

    // otherwise the preconditions for SSE codes are satisfied. 
    if(p==8)
    {
	// specific for p=8
      __DISPATCHER_MSG("[FGG INT SSE] P=8\n");
      SE_FKG_int_split_SSE_P8(phi, work, params);
    } 
    else if(p==16)
    {
    	// specific for p=16
    	__DISPATCHER_MSG("[FKG INT SSE] P=16\n");
	SE_FKG_int_split_SSE_P16(phi, work, params);
    }
    else if(p==6)
    {
    	// specific for p=6
    	__DISPATCHER_MSG("[FKG INT SSE] P=6\n");
    	SE_FKG_int_split_SSE_P6(phi, work, params);
    }    
    else
    {
      // vanilla SSE code (any even p)
      __DISPATCHER_MSG("[FKG INT SSE] Vanilla\n");
      SE_FKG_int_split_SSE(phi, work, params);
    }
}

void SE_FKG_int_split_AVX_dispatch(double* restrict phi,
                                   const SE_FGG_work* work,
                                   const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#ifdef AVX_FMA
    __DISPATCHER_MSG("[FKG INT AVX-FMA] ");
#else
    __DISPATCHER_MSG("[FGG INT AVX] ");
#endif
#if 0
    // THIS BYPASSES THE FAST AVX KERNELS.
    __DISPATCHER_MSG("AVX Disabled\n");
    SE_FKG_int_split(phi, work, params);
    return;
#endif
    
    // if either P or increments are not divisible by 4, fall back to SSE
    if( isnot_div_by_4(incri) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
        __DISPATCHER_MSG("AVX Abort (PARAMS)\n");
        SE_FKG_int_split_SSE_dispatch(phi, work, params);
        return;
    }

    // otherwise the preconditions for AVX codes are satisfied. 
    if(p==16)
    {
        // specific for p=16
        __DISPATCHER_MSG(" P=16\n");
        SE_FKG_int_split_AVX_P16(phi, work, params);
    }
    else if(p==8)
    {
        // specific for p=8
        __DISPATCHER_MSG("P=8\n");
        SE_FKG_int_split_AVX_P8(phi, work, params);
    }
    else if(p==4)
    {
        // specific for p=4
        __DISPATCHER_MSG("P=4\n");
        SE_FKG_int_split_AVX_P4(phi, work, params);
    }
    else
    {
      __DISPATCHER_MSG("AVX Abort\n");
      SE_FKG_int_split_SSE_dispatch(phi, work, params);
    }
}

void SE_FKG_grid(SE_FGG_work* work, const SE_state* st,
		 const SE_FGG_params* params)
{
    // vectors for FGG expansions
    double zx0[P_MAX] MEM_ALIGNED;
    double zy0[P_MAX] MEM_ALIGNED;
    double zz0[P_MAX] MEM_ALIGNED;

    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const int p = params->P;

    double cij0;
    double xn[3];
    double qn;
    int idx0, i,j,k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++) {
      // compute index and expansion vectors
      xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
      qn = st->q[n];
      
      idx0 = kaiser_expansion_3p(xn, qn, params, zx0, zy0, zz0);

      for(i = 0; i<p; i++) {
	  for(j = 0; j<p; j++) {
	    cij0 = zx0[i]*zy0[j];
	    for(k = 0; k<p; k++) {
	      H[idx0] += zz0[k]*cij0;
	      idx0++;
	    }
	    idx0 += incrj;
	  }
	  idx0 += incri;
      }
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split(SE_FGG_work* work, const SE_state* st, 
		       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double cij0;
    double qn;
    int idx0, idxzz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	qn = st->q[n];
	// inline vanilla loop
	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij0 = qn*zx[p*n+i]*zy[p*n+j];
		idxzz=p*n;
		for(k = 0; k<p; k++)
		{
		    H[idx0] += zz[idxzz]*cij0;
		    idx0++; idxzz++;
		}
		idx0 += incrj; 
	    }
	    idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_SSE(SE_FGG_work* work, const SE_state* st, 
			   const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rC;
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
      {
	idx0 = work->idx[n];

	qn = st->q[n];
	for(i = 0; i<p; i++)
	  {
	    for(j = 0; j<p; j++)
	      {
		rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		idx_zz=p*n;
		for(k = 0; k<p; k+=2)
		  {
		    rH0  = _mm_loadu_pd( H+idx0 );
		    rZZ0 = _mm_load_pd( zz + idx_zz );
#ifdef AVX_FMA
		    rH0  = _mm_fmadd_pd(rC,rZZ0,rH0);		    
#else
		    rZZ0 = _mm_mul_pd(rZZ0,rC);
		    rH0  = _mm_add_pd(rH0,rZZ0);		    
#endif
		    _mm_storeu_pd( H+idx0, rH0 );

		    idx0  +=2;
		    idx_zz+=2;
		  }
		idx0 += incrj;
	      }
	    idx0 += incri;
	  }
      }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_SSE_P6(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    double qn;
    int idx0, idx_zz, i, j;
    const int incrj = params->npdims[2]-6; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-6);// outer increment

    __m128d rH0, rZZ0, rC;
    __m128d rH1, rZZ1;
    __m128d rH2, rZZ2;
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
      {
	idx0 = work->idx[n];
	int sixn = 6*n;
	qn = st->q[n];
	for(i = 0; i<6; i++)
	  {
	    double tmp = qn*zx[sixn+i];
	    for(j = 0; j<6; j++)
	      {
		rC = _mm_set1_pd( tmp*zy[sixn+j] );
		idx_zz = sixn;

		rH0  = _mm_loadu_pd( H+idx0 );
		rH1  = _mm_loadu_pd( H+idx0 + 2);
		rH2  = _mm_loadu_pd( H+idx0 + 4);
		
		rZZ0 = _mm_load_pd( zz + idx_zz );
		rZZ1 = _mm_load_pd( zz + idx_zz + 2);
		rZZ2 = _mm_load_pd( zz + idx_zz + 4);
#ifdef AVX_FMA
		rH0 = _mm_fmadd_pd(rZZ0,rC,rH0);
		rH1 = _mm_fmadd_pd(rZZ1,rC,rH1);
		rH2 = _mm_fmadd_pd(rZZ2,rC,rH2);		
#else
		rH0 = _mm_add_pd(rH0,_mm_mul_pd(rZZ0,rC));
		rH1 = _mm_add_pd(rH1,_mm_mul_pd(rZZ1,rC));
		rH2 = _mm_add_pd(rH2,_mm_mul_pd(rZZ2,rC));
#endif		
		_mm_storeu_pd( H+idx0    , rH0 );
		_mm_storeu_pd( H+idx0 + 2, rH1 );
		_mm_storeu_pd( H+idx0 + 4, rH2 );

		idx0 += incrj + 6;
	      }
	    idx0 += incri;
	  }
      }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_SSE_u8(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rC;
    __m128d rH1, rZZ1;
    __m128d rH2, rZZ2;
    __m128d rH3, rZZ3;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];

	qn = st->q[n];
	_mm_prefetch( (void*) (H+idx0), _MM_HINT_T0);
	int np = n*p;
	for(i = 0; i<p; i++)
	  {
	    double tmp = qn*zx[np+i];
	    for(j = 0; j<p; j++)
	      {
		rC = _mm_set1_pd( tmp*zy[np+j] );
		idx_zz=np;
		    
		for(k = 0; k<p; k+=8)
		  {
		    rH0  = _mm_loadu_pd( H+idx0     );
		    rH1  = _mm_loadu_pd( H+idx0 + 2 );
		    rH2  = _mm_loadu_pd( H+idx0 + 4 );
		    rH3  = _mm_loadu_pd( H+idx0 + 6 );

		    rZZ0 = _mm_load_pd( zz + idx_zz     );
		    rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
		    rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
		    rZZ3 = _mm_load_pd( zz + idx_zz + 6 );
#ifdef AVX_FMA
		    rH0 = _mm_fmadd_pd(rZZ0,rC,rH0);
		    rH1 = _mm_fmadd_pd(rZZ1,rC,rH1);
		    rH2 = _mm_fmadd_pd(rZZ2,rC,rH2);
		    rH3 = _mm_fmadd_pd(rZZ3,rC,rH3);		    
#else		    
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(rZZ0,rC));
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(rZZ1,rC));
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(rZZ2,rC));
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(rZZ3,rC));
#endif
		    _mm_storeu_pd( H+idx0    , rH0 );
		    _mm_storeu_pd( H+idx0 + 2, rH1 );
		    _mm_storeu_pd( H+idx0 + 4, rH2 );
		    _mm_storeu_pd( H+idx0 + 6, rH3 );

		    idx0  +=8;
		    idx_zz+=8;
		  }
		idx0 += incrj;
	      }
	    idx0 += incri;
	  }
    }
}

void SE_FKG_grid_split_SSE_P16(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, i, j;
    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

    __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
    __m128d rH0, rH1, rH2, rH3;
    __m128d rC;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
	idx = work->idx[n];

	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	qn = st->q[n];

	int n16 = n*16;
	
        rZZ0 = _mm_load_pd(zz + n16     );
        rZZ1 = _mm_load_pd(zz + n16 + 2 );
        rZZ2 = _mm_load_pd(zz + n16 + 4 );
        rZZ3 = _mm_load_pd(zz + n16 + 6 );
        rZZ4 = _mm_load_pd(zz + n16 + 8 );
        rZZ5 = _mm_load_pd(zz + n16 + 10);
        rZZ6 = _mm_load_pd(zz + n16 + 12);
        rZZ7 = _mm_load_pd(zz + n16 + 14);

	for(i = 0; i<16; i++)
	    {
	      double tmp = qn*zx[n16+i];
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( tmp*zy[n16+j] );

                    /* 0 - 3 */ 
                    rH0  = _mm_loadu_pd( H+idx    );
                    rH1  = _mm_loadu_pd( H+idx + 2);
                    rH2  = _mm_loadu_pd( H+idx + 4);
                    rH3  = _mm_loadu_pd( H+idx + 6);

		    // if zs does not have 16-byte alignment, this will core.
		    // PLATFORM AND COMPILER DEPENDENT (FIXME)
#ifdef AVX_FMA
		    rH0 = _mm_add_pd(rZZ0,rC,rH0);
		    rH1 = _mm_add_pd(rZZ1,rC,rH1);
		    rH2 = _mm_add_pd(rZZ2,rC,rH2);
		    rH3 = _mm_add_pd(rZZ3,rC,rH3);
#else		    
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(rZZ0,rC));
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(rZZ1,rC));
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(rZZ2,rC));
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(rZZ3,rC));
#endif	
		    _mm_storeu_pd(H + idx    , rH0);
		    _mm_storeu_pd(H + idx + 2, rH1);
		    _mm_storeu_pd(H + idx + 4, rH2);
		    _mm_storeu_pd(H + idx + 6, rH3);

                    /* 4 - 7*/ 
		    rH0  = _mm_loadu_pd( H+idx + 8 );
                    rH1  = _mm_loadu_pd( H+idx + 10);
                    rH2  = _mm_loadu_pd( H+idx + 12);
                    rH3  = _mm_loadu_pd( H+idx + 14);

#ifdef AVX_FMA
		    rH0 = _mm_add_pd(rZZ4,rC,rH0);
		    rH1 = _mm_add_pd(rZZ5,rC,rH1);
		    rH2 = _mm_add_pd(rZZ6,rC,rH2);
		    rH3 = _mm_add_pd(rZZ7,rC,rH3);
#else
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(rZZ4,rC));
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(rZZ5,rC));
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(rZZ6,rC));
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(rZZ7,rC));
#endif	

		    _mm_storeu_pd(H + idx + 8 , rH0);
		    _mm_storeu_pd(H + idx + 10, rH1);
		    _mm_storeu_pd(H + idx + 12, rH2);
		    _mm_storeu_pd(H + idx + 14, rH3);

		    idx += incrj + 16;
		}
		idx += incri;
	    }
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_AVX_P4(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, i, j;
    const int incrj = params->npdims[2]-4; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-4);// outer increment

    __m256d rZZ0;
    __m256d rH0;
    __m256d rC;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
      {
	idx = work->idx[n];
	qn = st->q[n];
	int n4 = n*4;
        rZZ0 = _mm256_load_pd(zz + n4);

	for(i = 0; i<4; i++)
	  {
	    double tmp = qn*zx[n4+i];
	    for(j = 0; j<4; j++)
	      {
		rC = _mm256_set1_pd( tmp*zy[n4+j] );

		rH0  = _mm256_loadu_pd( H+idx     );

#ifdef AVX_FMA
		rH0 = _mm256_fmadd_pd(rZZ0,rC, rH0);
#else
		rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
#endif
		_mm256_storeu_pd(H + idx, rH0);

		idx += incrj + 4;
	      }
	    idx += incri;
	  }
      }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_AVX_P8(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, i, j;
    const int incrj = params->npdims[2]-8; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-8);// outer increment

    __m256d rZZ0, rZZ1; 
    __m256d rH0, rH1;
    __m256d rC;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
      {
	idx = work->idx[n];

	qn = st->q[n];
        rZZ0 = _mm256_load_pd(zz + n*8     );
        rZZ1 = _mm256_load_pd(zz + n*8 + 4 );

	for(i = 0; i<8; i++)
	  {
	    for(j = 0; j<8; j++)
	      {
		rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

		rH0  = _mm256_loadu_pd( H+idx     );
		rH1  = _mm256_loadu_pd( H+idx +  4);

#ifdef AVX_FMA
		rH0 = _mm256_fmadd_pd(rZZ0,rC, rH0);
		rH1 = _mm256_fmadd_pd(rZZ1,rC, rH1);
#else
		rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
		rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(rZZ1,rC));
#endif
		_mm256_storeu_pd(H + idx,      rH0);
		_mm256_storeu_pd(H + idx + 4,  rH1);

		idx += incrj + 8;
	      }
	    idx += incri;
	  }
      }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_AVX_P16(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;


    double qn;
    int idx, i, j;

    __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m256d rH0, rH1, rH2, rH3;
    __m256d rC;

    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
      {
	idx = work->idx[n];

	qn = st->q[n];
	int n16 = n*16;
	
        rZZ0 = _mm256_load_pd(zz + n16     );
        rZZ1 = _mm256_load_pd(zz + n16 + 4 );
        rZZ2 = _mm256_load_pd(zz + n16 + 8 );
        rZZ3 = _mm256_load_pd(zz + n16 + 12);

	
	for(i = 0 ; i<16; i++)
	  {
	    double tmp = qn*zx[n16+i];
	    for(j = 0; j<16; j++)
	      {
		rC = _mm256_set1_pd( tmp*zy[n16+j] );

		rH0  = _mm256_loadu_pd( H+idx     );
		rH1  = _mm256_loadu_pd( H+idx +  4);
		rH2  = _mm256_loadu_pd( H+idx +  8);
		rH3  = _mm256_loadu_pd( H+idx + 12);

#ifdef AVX_FMA
		rH0 = _mm256_fmadd_pd(rC,rZZ0,rH0);
		rH1 = _mm256_fmadd_pd(rC,rZZ1,rH1);
		rH2 = _mm256_fmadd_pd(rC,rZZ2,rH2);
		rH3 = _mm256_fmadd_pd(rC,rZZ3,rH3);
#else
		rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
		rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(rZZ1,rC));
		rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(rZZ2,rC));
		rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(rZZ3,rC));
#endif
		_mm256_storeu_pd(H + idx,      rH0);
		_mm256_storeu_pd(H + idx + 4,  rH1);
		_mm256_storeu_pd(H + idx + 8,  rH2);
		_mm256_storeu_pd(H + idx + 12, rH3);

		idx += incrj + 16;
	      }
	    idx += incri;
	  }
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_grid_split_SSE_dispatch(SE_FGG_work* work, const SE_state* st, 
				    const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment
  
#if 0
  __DISPATCHER_MSG("[FKG GRID SSE] SSE Disabled\n");
  SE_FKG_grid_split(work, st, params);
  return;
#endif
  
  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
      __DISPATCHER_MSG("[FKG GRID SSE] SSE Abort (PARAMS)\n");
      SE_FKG_grid_split(work, st, params);
      return;
    }
  
#if 0
  if( ( (unsigned long) work->zx)%16 != 0 || 
      ( (unsigned long) work->zy)%16 != 0 ||
      ( (unsigned long) work->zz)%16 != 0 )
    {
      __DISPATCHER_MSG("[FKG GRID SSE] SSE Abort (DATA)\n");
      SE_FKG_grid_split(work, st, params);
      return;
    }
#endif
  
  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==16)
    {
      // specific for p=16
      __DISPATCHER_MSG("[FKG GRID SSE] P=16\n");
      SE_FKG_grid_split_SSE_P16(work, st, params);
    }
  else if(p%8==0)
    {
      // specific for p divisible by 8
      __DISPATCHER_MSG("[FKG GRID SSE] P unroll 8\n");
      SE_FKG_grid_split_SSE_u8(work, st, params);
    }
  else if (p==6)
    {
      __DISPATCHER_MSG("[FKG GRID SSE] P=6\n");
      SE_FKG_grid_split_SSE_P6(work, st, params);
    }
  else
    {
      // vanilla SSE code (any even p)
      __DISPATCHER_MSG("[FKG GRID SSE] Vanilla\n");
      SE_FKG_grid_split_SSE(work, st, params);
    }
}

void SE_FKG_grid_split_AVX_dispatch(SE_FGG_work* work, const SE_state* st, 
				    const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#ifdef AVX_FMA
    __DISPATCHER_MSG("[FKG GRID AVX-FMA");
#else
    __DISPATCHER_MSG("[FKG GRID AVX]");
#endif

#if 0
    __DISPATCHER_MSG("AVX Disabled\n");
    SE_FKG_grid_split(work, st, params);
    return;
#endif

    // if either P or increments are not divisible by 4, fall back to SSE
    if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
	__DISPATCHER_MSG("AVX Abort (PARAMS)\n");
	SE_FKG_grid_split_SSE_dispatch(work, st, params);
	return;
    }

#if 0
    if(	( (unsigned long) work->zx)%32 != 0 || 
	( (unsigned long) work->zy)%32 != 0 ||
	( (unsigned long) work->zz)%32 != 0 )
    {
	__DISPATCHER_MSG("AVX Abort (DATA)\n");
	SE_FKG_grid_split(work, st, params);
	return;
    }
#endif

    // otherwise the preconditions for AVX codes are satisfied. 
    if(p==16)
    {
    	// specific for p=16
    	__DISPATCHER_MSG(" P=16\n");
    	SE_FKG_grid_split_AVX_P16(work, st, params);
    }
    else if(p==8)
    {
      // specific for p=8
      __DISPATCHER_MSG(" P=8\n");
      SE_FKG_grid_split_AVX_P8(work, st, params);
    }
    else if(p==4)
    {
      // specific for p=4
      __DISPATCHER_MSG(" P=4\n");
      SE_FKG_grid_split_AVX_P4(work, st, params);
    }      
    else
    {
      __DISPATCHER_MSG("AVX Abort\n");
      SE_FKG_grid_split_SSE_dispatch(work, st, params);
    }
}
