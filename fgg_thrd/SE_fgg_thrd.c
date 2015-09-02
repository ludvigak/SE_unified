#include "SE_fgg.h"
#include "SE_fgg_thrd.h"

#include "time.h"
#include "sys/time.h"
#include "unistd.h"

#include "fgg_thrd.h"
#include "fgg_thrd.c"

#define DELTA(tic,toc) ((toc.tv_sec  - tic.tv_sec) * 1000000u + toc.tv_usec - tic.tv_usec) / 1.e6

// Ludvig af Klinteberg, ludvigak@kth.se

// =============================================================================
//#define PI 3.141592653589793

// =============================================================================
// Internal routines ===========================================================
static 
int fgg_expansion_3p(const double x[3], const double q,
		     const SE_FGG_params* params,
		     double z2_0[P_MAX], 
		     double z2_1[P_MAX], 
		     double z2_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    
    double t0[3];
    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p))
    {
	for(int j=0; j<3; j++)
	{
	    idx = (int) round(x[j]/h);
	    idx_from[j] = idx - p_half;
	    t0[j] = x[j]-h*idx;
	}
    }
    else
    {
	for(int j=0; j<3; j++)
	{
	    idx = (int) floor(x[j]/h);
	    idx_from[j] = idx - (p_half-1);
	    t0[j] = x[j]-h*idx;
	}
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
    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;
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

// =============================================================================


// GRID_THRD ########################


// ##########################

void SE_FGG_grid_thrd(SE_FGG_work* work, const SE_state* st, 
		 const SE_FGG_params* params)
{
    // vectors for FGG expansions
    double zx0[P_MAX] MEM_ALIGNED;
    double zy0[P_MAX] MEM_ALIGNED;
    double zz0[P_MAX] MEM_ALIGNED;

    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const int p = params->P;
    
    double cij0;
    double xn[3];
    double qn;
    int idx0, zidx, i,j,k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    // ### GRID_THRD setup ###
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);
    // #######################

    // Loop over points
    for(int n=0; n<N; n++)
    {
	// compute index and expansion vectors
	xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
	qn = st->q[n];       
	idx0 = __FGG_EXPA(xn, qn, params, zx0, zy0, zz0);
	// ### GRID_THRD slice ###
	int i_end;	
	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &zidx);
	if (grid_thrd_ws.skip) 
	    continue;
	// #######################
	for(; i < i_end; i++)
	{	    
	    for(j = 0; j<p; j++)
	    {
		cij0 = zx0[i]*zy0[j];
		for(k = 0; k<p; k++)
		{
		    H[idx0] += zs[zidx]*zz0[k]*cij0;
		    idx0++; 
		    zidx++;
		}
		idx0 += incrj; 
	    }
	    idx0 += incri; 
	}
    }
}
