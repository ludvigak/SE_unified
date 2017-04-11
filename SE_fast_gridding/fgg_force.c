#ifdef FORCE
#ifdef THREE_PERIODIC
static 
int fgg_expansion_3p_force(const double x[3], const double q,
			   const SE_FGG_params* params,
			   double z2_0[P_MAX], 
			   double z2_1[P_MAX], 
			   double z2_2[P_MAX],
			   double zf_0[P_MAX],
			   double zf_1[P_MAX],
			   double zf_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    
    double t0[3];
    int idx;
    int idx_from[3],p_from;

    // compute index range and centering
    fgg_offset_3p(x, params, t0, idx_from);

    // compute third factor 
    double z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

    // compute second factor by induction
    double z_base0 = exp(2*c*h*t0[0]);
    double z_base1 = exp(2*c*h*t0[1]);
    double z_base2 = exp(2*c*h*t0[2]);

    double z0, z1, z2;
    if(is_odd(p)) 
    {
	z0 = pow(z_base0,-p_half);
	z1 = pow(z_base1,-p_half);
	z2 = pow(z_base2,-p_half);
	p_from = -p_half;
    }	
    else
    {
    	z0 = pow(z_base0,-p_half+1);
    	z1 = pow(z_base1,-p_half+1);
    	z2 = pow(z_base2,-p_half+1);
	p_from = -p_half+1;
    }

    z2_0[0] = z0;
    z2_1[0] = z1;
    z2_2[0] = z2;
    
    // extra terms multiplied to calculate forces
    zf_0[0] = -c*(t0[0]-p_from*h);
    zf_1[0] = -c*(t0[1]-p_from*h);
    zf_2[0] = -c*(t0[2]-p_from*h);

    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;

	zf_0[i] = -c*(t0[0]-(p_from+i)*h);
	zf_1[i] = -c*(t0[1]-(p_from+i)*h);
	zf_2[i] = -c*(t0[2]-(p_from+i)*h);
    }
 
    // save some flops by multiplying one vector with z3 factor
    for(int i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}
#elif ONE_PERIODIC
// -----------------------------------------------------------------------------
static 
int fgg_expansion_1p_force(const double x[3], const double q,
			   const SE_FGG_params* params,
			   double z2_0[P_MAX], 
			   double z2_1[P_MAX], 
			   double z2_2[P_MAX],
                           double zf_0[P_MAX],
                           double zf_1[P_MAX],
                           double zf_2[P_MAX])
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    const double a=params->a;
    const double b=params->b;

    double t0[3];
    int idx;
    int idx_from[3], p_from;

    // compute index range and centering
    if(is_odd(p))
    {
	idx = (int) round((x[0]-(a+h/2))/h);
	idx_from[0] = idx - p_half;
	t0[0] = x[0] - (idx*h + (a+h/2));

	idx = (int) round((x[1]-(b+h/2))/h);
	idx_from[1] = idx - p_half;
	t0[1] = x[1] - (idx*h + (b+h/2));

	idx = (int) round(x[2]/h);
	idx_from[2] = idx - p_half;
	t0[2] = x[2]-h*idx;

	p_from = -p_half;
    }
    else
    {
	idx = (int) floor((x[0]-(a+h/2))/h);
	idx_from[0] = idx - (p_half-1);
	t0[0] = x[0] - (idx*h + (a+h/2));

	idx = (int) floor((x[1]-(b+h/2))/h);
	idx_from[1] = idx - (p_half-1);
	t0[1] = x[1] - (idx*h + (b+h/2));

	idx = (int) floor(x[2]/h);
	idx_from[2] = idx - (p_half-1);
	t0[2] = x[2]-h*idx;

	p_from = -p_half + 1;
    }

    // compute third factor 
    double z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

    // compute second factor by induction
    double z_base0 = exp(2*c*h*t0[0]);
    double z_base1 = exp(2*c*h*t0[1]);
    double z_base2 = exp(2*c*h*t0[2]);

    double z0, z1, z2;
    if(is_odd(p))
    {
	z0 = pow(z_base0,-p_half);
	z1 = pow(z_base1,-p_half);
	z2 = pow(z_base2,-p_half);
    }	
    else
    {
    	z0 = pow(z_base0,-p_half+1);
    	z1 = pow(z_base1,-p_half+1);
    	z2 = pow(z_base2,-p_half+1);
    }

    z2_0[0] = z0;
    z2_1[0] = z1;
    z2_2[0] = z2;

    // extra terms multiplied to calculate forces
    zf_0[0] = -2*c*(t0[0]-p_from*h);
    zf_1[0] = -2*c*(t0[1]-p_from*h);
    zf_2[0] = -2*c*(t0[2]-p_from*h);

    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;

        zf_0[i] = -2*c*(t0[0]-(p_from+i)*h);
        zf_1[i] = -2*c*(t0[1]-(p_from+i)*h);
        zf_2[i] = -2*c*(t0[2]-(p_from+i)*h);
    }

    // save some flops by multiplying one vector with z3 factor
    for(int i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0], 
		       idx_from[1], 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}
#endif // end Periodicty

// -----------------------------------------------------------------------------
void SE_FGG_expand_all_force(SE_FGG_work* work, 
				 const SE_state* st, 
				 const SE_FGG_params* params)
{
    double xn[3] MEM_ALIGNED;
    
    const int N = params->N;
    const int P = params->P;

    for(int n=0; n<N; n++)
    {
	// compute index and expansion vectors
	xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
	
	*(work->idx+n) = __FGG_EXPA_FORCE(xn,1,params, 
					  work->zx+n*P, 
					  work->zy+n*P, 
					  work->zz+n*P,
					  work->zfx+n*P,
					  work->zfy+n*P,
					  work->zfz+n*P); 
    }
}

// -----------------------------------------------------------------------------
// vanilla grid gather to calculate forces
void SE_FGG_int_force(double* restrict force,  
		      const SE_FGG_work* work, 
		      SE_state* st, 
		      const SE_FGG_params* params)
{
    double z2_0[P_MAX] MEM_ALIGNED;
    double z2_1[P_MAX] MEM_ALIGNED;
    double z2_2[P_MAX] MEM_ALIGNED;

    // to alculate forces
    double zf_0[P_MAX] MEM_ALIGNED;
    double zf_1[P_MAX] MEM_ALIGNED;
    double zf_2[P_MAX] MEM_ALIGNED;

    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    double xm[3];
    int i,j,k,idx, zidx;
    double force_m[3], cij,Hzc;
#ifdef CALC_ENERGY
    double phi_m;
#endif
    
    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
      xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N]; 

      idx = __FGG_EXPA_FORCE(xm, 1, params, z2_0, z2_1, z2_2,zf_0,zf_1,zf_2);
      
      force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;

#ifdef CALC_ENERGY
      phi_m = 0;
#endif
      zidx = 0;
      
      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij = z2_0[i]*z2_1[j];
	      for(k = 0; k<p; k++)
		{
		  Hzc         = H[idx]*zs[zidx]*z2_2[k]*cij;
#ifdef CALC_ENERGY
		  phi_m      += Hzc;
#endif
		  force_m[0] += Hzc*zf_0[i];
		  force_m[1] += Hzc*zf_1[j];
		  force_m[2] += Hzc*zf_2[k];	 

		  idx++; zidx++;
		}
	      idx += incrj;
	    }
	  idx += incri;
	}
      force[m    ] = h3*force_m[0];
      force[m+  N] = h3*force_m[1];
      force[m+2*N] = h3*force_m[2];
#ifdef CALC_ENERGY
      st->phi[m]    = h3*phi_m;
#endif
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_dispatch_force(double* restrict force,  
					 SE_state *st,
					 const SE_FGG_work* work, 
					 const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    // THIS BYPASSES THE FAST SSE KERNELS.
    __DISPATCHER_MSG("[FGG INT SSE] SSE Disabled\n");
    SE_FGG_int_split_force(force, st, work, params);
    return;
#endif

    // if P is odd, or if either increment is odd, fall back on vanilla
    if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
	__DISPATCHER_MSG("[FGG INT SSE] SSE Abort (PARAMS)\n");
	SE_FGG_int_split_force(force, st, work, params);
	return;
    }
    
    // otherwise the preconditions for SSE codes are satisfied. 
    if(p==8)
    {
	// specific for p=8
	__DISPATCHER_MSG("[FGG INT SSE] P=8\n");
	SE_FGG_int_split_SSE_P8_force(force, st, work, params);
    }
    else if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG INT SSE] P=16\n");
	SE_FGG_int_split_SSE_P16_force(force, st, work, params); 
    }
    else if(p%8==0)
    {
	// for p divisible by 8
	__DISPATCHER_MSG("[FGG INT SSE] P unroll 8\n");
	SE_FGG_int_split_SSE_u8_force(force, st, work, params); 
    }
    else
    {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG INT SSE] Vanilla\n");
	SE_FGG_int_split_SSE_force(force, st, work, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_force(double* restrict force,  
			    SE_state* st,
			    const SE_FGG_work* work, 
			    const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H   = work->H;
    const double* restrict zs  = work->zs;
    const double* restrict zx  = work->zx;
    const double* restrict zy  = work->zy;
    const double* restrict zz  = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;


    const int    p = params->P;
    const int    N = params->N;
    const double h = params->h;
    const double h3= h*h*h;

    int i,j,k,m,idx,idx_zs,idx_zz;
    double force_m[3], cij, Hzc;
#ifdef CALC_ENERGY
    double phi_m;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for private(m)// work-share over OpenMP threads here
#endif
    for(m=0; m<N; m++)
    {
	idx = work->idx[m];

	force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
        phi_m = 0;
#endif
	idx_zs = 0;

	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij = zx[m*p+i]*zy[m*p+j];
		idx_zz=m*p;
		for(k = 0; k<p; k++)
		  {
		    Hzc         = H[idx]*zs[idx_zs]*zz[idx_zz]*cij;
#ifdef CALC_ENERGY
		    phi_m      += Hzc;
#endif
		    force_m[0] += Hzc*zfx[m*p+i];
		    force_m[1] += Hzc*zfy[m*p+j];
		    force_m[2] += Hzc*zfz[m*p+k];
		    
		    idx++; idx_zs++; idx_zz++;
		}
		idx += incrj;
	    }
	    idx += incri;
	}
	force[m    ] = h3*force_m[0];
	force[m+  N] = h3*force_m[1];
	force[m+2*N] = h3*force_m[2];
#ifdef CALC_ENERGY
	st->phi[m]   = h3*phi_m;
#endif
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_force(double* restrict force,
				SE_state* st,
				const SE_FGG_work* work, 
				const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;


    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,k,m,idx,idx_zs,idx_zz;
    double sx[2] MEM_ALIGNED;
    double sy[2] MEM_ALIGNED;
    double sz[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rZS0,rZFZ0;
    __m128d rC, rCX, rCY;
    __m128d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
    double s[2]  MEM_ALIGNED;
    __m128d rP;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

    for(m=0; m<N; m++)
    {

	idx = work->idx[m];	
	idx_zs = 0;
	rFX = _mm_setzero_pd();
	rFY = _mm_setzero_pd();
	rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
	rP = _mm_setzero_pd();
#endif

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		  double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm_set1_pd( tmp );
		    rCX = _mm_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=2)
		    {
		      rZFZ0= _mm_load_pd( zfz+ m*p+k );
		      rH0  = _mm_load_pd( H  + idx );
		      rZZ0 = _mm_load_pd( zz + idx_zz);
		      rZS0 = _mm_load_pd( zs + idx_zs);
		      rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		      rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		      rFZ = _mm_add_pd(rFZ,_mm_mul_pd(_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
#endif			
			idx+=2; 
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	  for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		  double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm_set1_pd( tmp );
		    rCX = _mm_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=2)
		    {
		      rZFZ0= _mm_load_pd( zfz + m*p+k );
		      rH0  = _mm_loadu_pd( H+idx );
		      rZZ0 = _mm_load_pd( zz + idx_zz);
		      rZS0 = _mm_load_pd( zs + idx_zs);
		      rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		      rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		      rFZ = _mm_add_pd(rFZ,_mm_mul_pd(_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
#endif			
			idx+=2; 
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }

	}
	_mm_store_pd(sx,rFX);
	_mm_store_pd(sy,rFY);
	_mm_store_pd(sz,rFZ);
	
	force[m    ] = h3*(sx[0]+sx[1]);
	force[m+  N] = h3*(sy[0]+sy[1]);
	force[m+2*N] = h3*(sz[0]+sz[1]);

#ifdef CALC_ENERGY
	_mm_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]);
#endif

    }
}


// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P8_force(double* restrict force,  
				   SE_state* st,
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    /* ASSUME P=8 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx,idx_zs;

    double sx[2] MEM_ALIGNED;
    double sy[2] MEM_ALIGNED;
    double sz[2] MEM_ALIGNED;

    // hold entire zz vector
    __m128d rZZ0, rZZ1, rZZ2, rZZ3;
    __m128d rC, rCX, rCY;
    __m128d rH0, rH1, rH2, rH3; 
    __m128d rZS0, rZS1, rZS2, rZS3;
    __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
    __m128d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
    double s[2]  MEM_ALIGNED;
   __m128d rP;
#endif

    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);

    for(int m=0; m<N; m++)
    {

	idx = work->idx[m];
	idx_zs = 0;
	rFX = _mm_setzero_pd();
	rFY = _mm_setzero_pd();
	rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm_setzero_pd();
#endif


	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*8     );
	rZZ1 = _mm_load_pd(zz + m*8 + 2 );
	rZZ2 = _mm_load_pd(zz + m*8 + 4 );
	rZZ3 = _mm_load_pd(zz + m*8 + 6 );

	/* hoist load of ZFZ vector */
	rZFZ0 = _mm_load_pd(zfz + m*8     );
	rZFZ1 = _mm_load_pd(zfz + m*8 + 2 );
	rZFZ2 = _mm_load_pd(zfz + m*8 + 4 );
	rZFZ3 = _mm_load_pd(zfz + m*8 + 6 );

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<8; i++)
	    {
	      for(j = 0; j<8; j++)
		{
		  double tmp = zx[m*8+i]*zy[m*8 + j];
		  rC  = _mm_set1_pd( tmp );
		  rCX = _mm_set1_pd( tmp * zfx[m*8 + i]);
		  rCY = _mm_set1_pd( tmp * zfy[m*8 + j]);
		  
		  rH0  = _mm_load_pd( H+idx    );
		  rH1  = _mm_load_pd( H+idx + 2);
		  rH2  = _mm_load_pd( H+idx + 4);
		  rH3  = _mm_load_pd( H+idx + 6);
		  
		  rZS0 = _mm_load_pd( zs + idx_zs    );
		  rZS1 = _mm_load_pd( zs + idx_zs + 2);
		  rZS2 = _mm_load_pd( zs + idx_zs + 4);
		  rZS3 = _mm_load_pd( zs + idx_zs + 6);
		  
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));		  

		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));
		  
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));
		  
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
#endif


		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	  for(i = 0; i<8; i++)
	    {
	      for(j = 0; j<8; j++)
		{
		  double tmp = zx[m*8+i]*zy[m*8 + j];
		  rC  = _mm_set1_pd( tmp );
		  rCX = _mm_set1_pd( tmp * zfx[m*8 + i]);
		  rCY = _mm_set1_pd( tmp * zfy[m*8 + j]);

		  rH0  = _mm_loadu_pd( H+idx    );
		  rH1  = _mm_loadu_pd( H+idx + 2);
		  rH2  = _mm_loadu_pd( H+idx + 4);
		  rH3  = _mm_loadu_pd( H+idx + 6);

		  rZS0 = _mm_load_pd( zs + idx_zs    );
		  rZS1 = _mm_load_pd( zs + idx_zs + 2);
		  rZS2 = _mm_load_pd( zs + idx_zs + 4);
		  rZS3 = _mm_load_pd( zs + idx_zs + 6);
		    
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));

		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));

		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
#endif

		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	_mm_store_pd(sx,rFX);
	_mm_store_pd(sy,rFY);
	_mm_store_pd(sz,rFZ);

	force[m    ] = h3*(sx[0]+sx[1]);
	force[m+  N] = h3*(sy[0]+sy[1]);
	force[m+2*N] = h3*(sz[0]+sz[1]);

#ifdef CALC_ENERGY
	_mm_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]);
#endif

    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P16_force(double* restrict force,  
			            SE_state* st,
			      const SE_FGG_work* work, 
			      const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,idx,idx_zs;
    double sx[2] MEM_ALIGNED;
    double sy[2] MEM_ALIGNED;
    double sz[2] MEM_ALIGNED;
    

    // hold entire zz vector
    __m128d rZZ0 , rZZ1 , rZZ2 , rZZ3 , rZZ4 , rZZ5 , rZZ6 , rZZ7; 
    __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3, rZFZ4, rZFZ5, rZFZ6, rZFZ7;
    __m128d rC, rCX, rCY, rFX, rFY, rFZ;
    __m128d rH0, rZS0;

#ifdef CALC_ENERGY
    double s[2]  MEM_ALIGNED;
    __m128d rP;
#endif

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);

    for(int m=0; m<N; m++)
      {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	rFX = _mm_setzero_pd();
	rFY = _mm_setzero_pd();
	rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm_setzero_pd();
#endif

	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*16     );
	rZZ1 = _mm_load_pd(zz + m*16 + 2 );
	rZZ2 = _mm_load_pd(zz + m*16 + 4 );
	rZZ3 = _mm_load_pd(zz + m*16 + 6 );
	rZZ4 = _mm_load_pd(zz + m*16 + 8 );
	rZZ5 = _mm_load_pd(zz + m*16 + 10);
	rZZ6 = _mm_load_pd(zz + m*16 + 12);
	rZZ7 = _mm_load_pd(zz + m*16 + 14);

	/* hoist load of ZFZ vector */
	rZFZ0 = _mm_load_pd(zfz + m*16     );
	rZFZ1 = _mm_load_pd(zfz + m*16 + 2 );
	rZFZ2 = _mm_load_pd(zfz + m*16 + 4 );
	rZFZ3 = _mm_load_pd(zfz + m*16 + 6 );
	rZFZ4 = _mm_load_pd(zfz + m*16 + 8 );
	rZFZ5 = _mm_load_pd(zfz + m*16 + 10);
	rZFZ6 = _mm_load_pd(zfz + m*16 + 12);
	rZFZ7 = _mm_load_pd(zfz + m*16 + 14);

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<16; i++)
	    {
	      for(j = 0; j<16; j++)
		{
		  double tmp = zx[m*16+i]*zy[m*16+j];
		  rC  = _mm_set1_pd( tmp );
		  rCX = _mm_set1_pd( tmp * zfx[m*16+i] );
		  rCY = _mm_set1_pd( tmp * zfy[m*16+j]);

		  /* 0 */ 
		  rH0  = _mm_load_pd( H+idx );
		  rZS0 = _mm_load_pd( zs + idx_zs);		    
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));

#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
#endif
		  /* 1 */ 
		  rH0  = _mm_load_pd( H+idx + 2);
		  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));
#endif

		  /* 2 */ 
		  rH0  = _mm_load_pd( H+idx + 4);
		  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));
#endif

		  /* 3 */ 
		  rH0  = _mm_load_pd( H+idx + 6);
		  rZS0 = _mm_load_pd( zs + idx_zs + 6);
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));
#endif

		  /* 4 */ 
		  rH0  = _mm_load_pd( H+idx + 8);
		  rZS0 = _mm_load_pd( zs + idx_zs + 8);
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_NERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));
#endif

		  /* 5 */ 
		  rH0  = _mm_load_pd( H+idx + 10);
		  rZS0 = _mm_load_pd( zs + idx_zs + 10);
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));
#endif

		  /* 6 */ 
		  rH0  = _mm_load_pd( H+idx + 12);
		  rZS0 = _mm_load_pd( zs + idx_zs + 12);
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));
#endif
		  /* 7 */ 
		  rH0  = _mm_load_pd( H+idx + 14);
		  rZS0 = _mm_load_pd( zs + idx_zs + 14);
		  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
		  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
		  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
		  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));
#endif

		  idx_zs +=16;
		  idx += incrj + 16;
		}
	      idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	  {
	    for(i = 0; i<16; i++)
	      {
		for(j = 0; j<16; j++)
		  {
		    double tmp = zx[m*16+i]*zy[m*16+j];
		    rC  = _mm_set1_pd( tmp );
                    rCX = _mm_set1_pd( tmp * zfx[m*16+i]);
                    rCY = _mm_set1_pd( tmp * zfy[m*16+j]);
		    
		    /* 0 */ 
		    rH0  = _mm_loadu_pd( H+idx );
		    rZS0 = _mm_load_pd( zs + idx_zs);		    
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
#endif

		    /* 1 */ 
		    rH0  = _mm_loadu_pd( H+idx + 2);
		    rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));
#endif

		    /* 2 */ 
		    rH0  = _mm_loadu_pd( H+idx + 4);
		    rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));
#endif

		    /* 3 */ 
		    rH0  = _mm_loadu_pd( H+idx + 6);
		    rZS0 = _mm_load_pd( zs + idx_zs + 6);
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));
#endif

		    /* 4 */ 
		    rH0  = _mm_loadu_pd( H+idx + 8);
		    rZS0 = _mm_load_pd( zs + idx_zs + 8);
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));
#endif

		    /* 5 */ 
		    rH0  = _mm_loadu_pd( H+idx + 10);
		    rZS0 = _mm_load_pd( zs + idx_zs + 10);
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));
#endif

		    /* 6 */ 
		    rH0  = _mm_loadu_pd( H+idx + 12);
		    rZS0 = _mm_load_pd( zs + idx_zs + 12);
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));
#endif

		    /* 7 */ 
		    rH0  = _mm_loadu_pd( H+idx + 14);
		    rZS0 = _mm_load_pd( zs + idx_zs + 14);
		    rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
		    rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
		    rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
		    rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));
#endif
		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}

	_mm_store_pd(sx,rFX);
	_mm_store_pd(sy,rFY);
	_mm_store_pd(sz,rFZ);

	force[m    ] = h3*(sx[0]+sx[1]);
	force[m+  N] = h3*(sy[0]+sy[1]);
	force[m+2*N] = h3*(sz[0]+sz[1]);

#ifdef CALC_ENERGY
	_mm_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]);
#endif
    }
}
// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_u8_force(double* restrict force,  
				   SE_state* st,
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;    

    int i,j,k,idx,idx_zs,idx_zz;
    double sx[2] MEM_ALIGNED;
    double sy[2] MEM_ALIGNED;
    double sz[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rZS0, rZFZ0;
    __m128d rH1, rZZ1, rZS1, rZFZ1;
    __m128d rH2, rZZ2, rZS2, rZFZ2;
    __m128d rH3, rZZ3, rZS3, rZFZ3;
    __m128d rFX, rFY, rFZ;
    __m128d  rC, rCX, rCY;

#ifdef CALC_ENERGY
    double s[2]  MEM_ALIGNED;
    __m128d rP;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

    for(int m=0; m<N; m++)
      {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	rFX = _mm_setzero_pd();
	rFY = _mm_setzero_pd();
	rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm_setzero_pd();
#endif

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		  {
		    double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm_set1_pd( tmp );
		    rCX = _mm_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		      {
			rH0  = _mm_load_pd( H+idx    );
			rH1  = _mm_load_pd( H+idx + 2);
			rH2  = _mm_load_pd( H+idx + 4);
			rH3  = _mm_load_pd( H+idx + 6);

			rZZ0 = _mm_load_pd( zz + idx_zz    );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2);
			rZZ2 = _mm_load_pd( zz + idx_zz + 4);
			rZZ3 = _mm_load_pd( zz + idx_zz + 6);

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rZFZ0 = _mm_load_pd(zfz+ idx_zz    );
			rZFZ1 = _mm_load_pd(zfz+ idx_zz + 2);
			rZFZ2 = _mm_load_pd(zfz+ idx_zz + 4);
			rZFZ3 = _mm_load_pd(zfz+ idx_zz + 6);

			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));

			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));
			
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
#endif

			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned load from H
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		  double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm_set1_pd( tmp );
		    rCX = _mm_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm_loadu_pd( H+idx    );
			rH1  = _mm_loadu_pd( H+idx + 2);
			rH2  = _mm_loadu_pd( H+idx + 4);
			rH3  = _mm_loadu_pd( H+idx + 6);

			rZZ0 = _mm_load_pd( zz + idx_zz    );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2);
			rZZ2 = _mm_load_pd( zz + idx_zz + 4);
			rZZ3 = _mm_load_pd( zz + idx_zz + 6);

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rZFZ0 = _mm_load_pd(zfz+ idx_zz    );
			rZFZ1 = _mm_load_pd(zfz+ idx_zz + 2);
			rZFZ2 = _mm_load_pd(zfz+ idx_zz + 4);
			rZFZ3 = _mm_load_pd(zfz+ idx_zz + 6);

			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
			rFX = _mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));

			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
			rFY = _mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));
			
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
			rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
#endif
			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}

	// done accumulating
	_mm_store_pd(sx,rFX);
	_mm_store_pd(sy,rFY);
	_mm_store_pd(sz,rFZ);
	
	force[m    ] = h3*(sx[0]+sx[1]);
	force[m+  N] = h3*(sy[0]+sy[1]);
	force[m+2*N] = h3*(sz[0]+sz[1]);

#ifdef CALC_ENERGY
	_mm_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]);
#endif
    }
}



// -----------------------------------------------------------------------------
#ifdef __AVX__
void SE_FGG_int_split_AVX_dispatch_force(double* restrict force,  
					 SE_state *st,
					 const SE_FGG_work* work, 
					 const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 1
    // THIS BYPASSES THE FAST AVX KERNELS.
    __DISPATCHER_MSG("[FGG INT AVX] AVX Disabled\n");
    SE_FGG_int_split_force(force, st, work, params);
    return;
#endif

    // if P, incri or increments are not divisible by 4, fall back on vanilla
    if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
	__DISPATCHER_MSG("[FGG INT AVX] AVX Abort (PARAMS)\n");
	SE_FGG_int_split_force(force, st, work, params);
	return;
    }
    
    // otherwise the preconditions for AVX codes are satisfied. 
    if(p==8)
    {
	// specific for p=8
	__DISPATCHER_MSG("[FGG INT AVX] P=8\n");
	SE_FGG_int_split_AVX_P8_force(force, st, work, params);
    }
    else if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG INT AVX] P=16\n");
	SE_FGG_int_split_AVX_P16_force(force, st, work, params); 
    }
    else if(p%8==0)
    {
	// for p divisible by 8
	__DISPATCHER_MSG("[FGG INT AVX] P unroll 8\n");
	SE_FGG_int_split_AVX_u8_force(force, st, work, params); 
    }
    else if(p%4==0)
    {
	// vanilla AVX code (p divisible by 4)
	__DISPATCHER_MSG("[FGG INT AVX] P unroll 4\n");
	SE_FGG_int_split_AVX_force(force, st, work, params);
    }
    else
    {
      // vanilla SSE code (any even p)
      __DISPATCHER_MSG("[FGG INT AVX] Vanilla\n");
      SE_FGG_int_split_SSE_force(force, st, work, params);
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_force(double* restrict force,
				SE_state* st,
				const SE_FGG_work* work, 
				const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;


    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3 = h*h*h;

    int i,j,k,m,idx,idx_zs,idx_zz;

    double sx[4] MEM_ALIGNED;
    double sy[4] MEM_ALIGNED;
    double sz[4] MEM_ALIGNED;

    __m256d rH0, rZZ0, rZS0,rZFZ0;
    __m256d rC, rCX, rCY;
    __m256d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
    double s[4]  MEM_ALIGNED;
    __m256d rP;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

    for(m=0; m<N; m++)
      {
	idx = work->idx[m];	
	idx_zs = 0;
	rFX = _mm256_setzero_pd();
	rFY = _mm256_setzero_pd();
	rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
	rP = _mm256_setzero_pd();
#endif

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		  {
		    double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=4)
		    {
		      rZFZ0= _mm256_load_pd( zfz+ m*p+k );
		      rH0  = _mm256_load_pd( H  + idx );
		      rZZ0 = _mm256_load_pd( zz + idx_zz);
		      rZS0 = _mm256_load_pd( zs + idx_zs);
		      rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		      rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		      rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
		      rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
#endif			
			idx+=4; 
			idx_zs+=4; 
			idx_zz+=4;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		  {
		    double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=4)
		    {
		      rZFZ0= _mm256_load_pd( zfz + m*p+k );
		      rH0  = _mm256_loadu_pd( H+idx );
		      rZZ0 = _mm256_load_pd( zz + idx_zz);
		      rZS0 = _mm256_load_pd( zs + idx_zs);
		      rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		      rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		      rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
		      rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
#endif			
			idx+=4; 
			idx_zs+=4; 
			idx_zz+=4;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }

	}
	_mm256_store_pd(sx,rFX);
	_mm256_store_pd(sy,rFY);
	_mm256_store_pd(sz,rFZ);
	force[m    ] = h3*(sx[0]+sx[1]+sx[2]+sx[3]);
	force[m+  N] = h3*(sy[0]+sy[1]+sy[2]+sy[3]);
	force[m+2*N] = h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
	_mm256_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
#endif

    }
}


// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_P16_force(double* restrict force,  
				   SE_state* st,
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3=h*h*h;
   
    int i,j,idx,idx_zs;

    double sx[4] MEM_ALIGNED;
    double sy[4] MEM_ALIGNED;
    double sz[4] MEM_ALIGNED;


    // hold entire zz vector
    __m256d rZZ0, rZZ1, rZZ2, rZZ3;
    __m256d rC, rCX, rCY;
    __m256d rH0, rH1, rH2, rH3; 
    __m256d rZS0, rZS1, rZS2, rZS3;
    __m256d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
    __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
    double s[4]  MEM_ALIGNED;
   __m256d rP;
#endif

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);


    for(int m=0; m<N; m++)
      {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0 );
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0 );
	rFX = _mm256_setzero_pd();
	rFY = _mm256_setzero_pd();
	rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm256_setzero_pd();
#endif


	/* hoist load of ZZ vector */
	rZZ0 = _mm256_load_pd(zz + m*16     );
	rZZ1 = _mm256_load_pd(zz + m*16 + 4 );
	rZZ2 = _mm256_load_pd(zz + m*16 + 8 );
	rZZ3 = _mm256_load_pd(zz + m*16 + 12);

	/* hoist load of ZFZ vector */
	rZFZ0 = _mm256_load_pd(zfz + m*16     );
	rZFZ1 = _mm256_load_pd(zfz + m*16 + 4 );
	rZFZ2 = _mm256_load_pd(zfz + m*16 + 8 );
	rZFZ3 = _mm256_load_pd(zfz + m*16 + 12);

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		  {
		    double tmp = zx[m*16+i]*zy[m*16 + j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*16 + i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*16 + j]);

		    rH0  = _mm256_load_pd( H+idx     );
		    rH1  = _mm256_load_pd( H+idx + 4 );
		    rH2  = _mm256_load_pd( H+idx + 8 );
		    rH3  = _mm256_load_pd( H+idx + 12);

		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		    rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
		    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCX),rZS2)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCX),rZS3)));

		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCY),rZS2)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCY),rZS3)));

		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3)));
#endif

		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 32-aligned, so use non-aligned loads
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		  {
		    double tmp = zx[m*16+i]*zy[m*16 + j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*16 + i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*16 + j]);

		    rH0  = _mm256_loadu_pd( H+idx     );
		    rH1  = _mm256_loadu_pd( H+idx + 4 );
		    rH2  = _mm256_loadu_pd( H+idx + 8 );
		    rH3  = _mm256_loadu_pd( H+idx + 12);

		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		    rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
		    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCX),rZS2)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCX),rZS3)));

		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCY),rZS2)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCY),rZS3)));

		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3)));
#endif

		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	_mm256_store_pd(sx,rFX);
	_mm256_store_pd(sy,rFY);
	_mm256_store_pd(sz,rFZ);

	force[m    ] = h3*(sx[0]+sx[1]+sx[2]+sx[3]);
	force[m+  N] = h3*(sy[0]+sy[1]+sy[2]+sy[3]);
	force[m+2*N] = h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
	_mm256_stream_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
#endif

    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_P8_force(double* restrict force,  
				   SE_state* st,
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    /* ASSUME P=8 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;
    const double h3=h*h*h;

    int i,j,idx,idx_zs;

    double sx[4] MEM_ALIGNED;
    double sy[4] MEM_ALIGNED;
    double sz[4] MEM_ALIGNED;

    // hold entire zz vector
    __m256d rZZ0, rZZ1;
    __m256d rC, rCX, rCY;
    __m256d rH0, rH1; 
    __m256d rZS0, rZS1;
    __m256d rZFZ0, rZFZ1;
    __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
    double s[4]  MEM_ALIGNED;
   __m256d rP;
#endif

    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);

    for(int m=0; m<N; m++)
      {
	idx = work->idx[m];
	idx_zs = 0;
	rFX = _mm256_setzero_pd();
	rFY = _mm256_setzero_pd();
	rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm256_setzero_pd();
#endif
	/* hoist load of ZZ vector */
	rZZ0 = _mm256_load_pd(zz + m*8     );
	rZZ1 = _mm256_load_pd(zz + m*8 + 4 );

	/* hoist load of ZFZ vector */
	rZFZ0 = _mm256_load_pd(zfz + m*8     );
	rZFZ1 = _mm256_load_pd(zfz + m*8 + 4 );

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		  {
		    double tmp = zx[m*8+i]*zy[m*8 + j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp  * zfx[m*8 + i]);
		    rCY = _mm256_set1_pd( tmp  * zfy[m*8 + j]);

		    rH0  = _mm256_load_pd( H+idx     );
		    rH1  = _mm256_load_pd( H+idx + 4 );
		 
		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		 		    
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		 

		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		 

		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		 

#ifdef CALC_ENERGY
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
		 
#endif
		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 32-aligned, so use non-aligned loads
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		  {
		    double tmp = zx[m*8+i]*zy[m*8 + j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp  * zfx[m*8 + i]);
		    rCY = _mm256_set1_pd( tmp  * zfy[m*8 + j]);

		    rH0  = _mm256_loadu_pd( H+idx     );
		    rH1  = _mm256_loadu_pd( H+idx + 4 );
		 
		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		 		    
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
		    rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		 
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
		    rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		 
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
		    rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		 

#ifdef CALC_ENERGY
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif

		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	_mm256_store_pd(sx,rFX);
	_mm256_store_pd(sy,rFY);
	_mm256_store_pd(sz,rFZ);

	force[m    ] = h3*(sx[0]+sx[1]+sx[2]+sx[3]);
	force[m+  N] = h3*(sy[0]+sy[1]+sy[2]+sy[3]);
	force[m+2*N] = h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
	_mm256_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
#endif

    }
}



// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_u8_force(double* restrict force,  
				   SE_state* st,
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;
    const double h3=h*h*h;

    int i,j,k,idx,idx_zs,idx_zz;

    double sx[4] MEM_ALIGNED;
    double sy[4] MEM_ALIGNED;
    double sz[4] MEM_ALIGNED;

    __m256d rH0, rZZ0, rZS0, rZFZ0;
    __m256d rH1, rZZ1, rZS1, rZFZ1;
    __m256d rFX, rFY, rFZ;
    __m256d  rC, rCX, rCY;

#ifdef CALC_ENERGY
    double s[4]  MEM_ALIGNED;
    __m256d rP;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

    for(int m=0; m<N; m++)
      {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);
	rFX = _mm256_setzero_pd();
	rFY = _mm256_setzero_pd();
	rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
	rP  = _mm256_setzero_pd();
#endif

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		  {
		    double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_load_pd( H+idx    );
			rH1  = _mm256_load_pd( H+idx + 4);

			rZZ0 = _mm256_load_pd( zz + idx_zz    );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);

			rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
			rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

			rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
			rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));

			rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
			rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
			
rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));

#ifdef CALC_ENERGY
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif

			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 32-aligned, so use non-aligned load from H
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		  {
		    double tmp = zx[m*p+i]*zy[m*p+j];
		    rC  = _mm256_set1_pd( tmp );
		    rCX = _mm256_set1_pd( tmp * zfx[m*p+i]);
		    rCY = _mm256_set1_pd( tmp * zfy[m*p+j]);

		    idx_zz=m*p;
		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_loadu_pd( H+idx    );
			rH1  = _mm256_loadu_pd( H+idx + 4);
		
			rZZ0 = _mm256_load_pd( zz + idx_zz    );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);

			rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
			rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

			rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
			rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));

			rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
			rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
			
			rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
			rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));

#ifdef CALC_ENERGY
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif

			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}

	// done accumulating
	_mm256_store_pd(sx,rFX);
	_mm256_store_pd(sy,rFY);
	_mm256_store_pd(sz,rFZ);

	force[m    ] = h3*(sx[0]+sx[1]+sx[2]+sx[3]);
	force[m+  N] = h3*(sy[0]+sy[1]+sy[2]+sy[3]);
	force[m+2*N] = h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
	_mm256_store_pd(s,rP);
	st->phi[m] = h3*(s[0]+s[1]+s[2]+s[3]);
#endif
    }
}
#endif // AVX

#endif  //end FORCE
