#include "SE_fgg.h"
#include "SE_fgg_thrd.h"

#include "time.h"
#include "sys/time.h"
#include "unistd.h"

#include "omp.h"

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
typedef struct {
    int p;    
    int block_begin;
    int block_end;
    int block_size;
    int skip;
} grid_thrd_ws_t;

static void inline grid_thrd_setup(grid_thrd_ws_t* ws, const int* npdims, int p)
{
#ifdef FGG_THRD
    // Determine which block current thread is working on
    int num_threads = omp_get_num_threads();
    int thread_num = omp_get_thread_num();
    int block_size = npdims[0] / num_threads;
    ws->block_begin = thread_num * block_size;
    ws->block_end = ws->block_begin+block_size;
    if (thread_num+1 == num_threads)
	ws->block_end = npdims[0];
    ws->block_size = npdims[1]*npdims[2];    
    ws->p = p;
#else
    ws->p = p;
#endif
}


static void inline grid_thrd_slice(grid_thrd_ws_t* ws, 
				   int* support_begin, 
				   int* first_slice, 
				   int* last_slice, 
				   int* zidx)
{
#ifdef FGG_THRD
    // Figure out which slices of current P^3 block that go into this block.
    *first_slice = 0;	
    *zidx = 0;
    int block_id = *support_begin / ws->block_size;
    // Skip point if above our chunk
    if (block_id >= ws->block_end)
    {
	ws->skip = true;
	return;
    }
    int diff = ws->block_begin - block_id;
    if (diff > 0)
    {
	// Skip point if too far below our chunk	    
	if (diff > ws->p)
	{
	    ws->skip = true;
	    return;
	}
	// Jump forward to first block inside our chunk
	*support_begin += diff*ws->block_size;
	*zidx += diff*ws->p*ws->p;
	block_id += diff;
	*first_slice += diff;
    }
    // Only iterate over blocks left inside our chunk, but max p of them
    *last_slice = *first_slice + (ws->block_end - block_id);
    if (*last_slice > ws->p)
	*last_slice = ws->p;
    ws->skip = false;
#else
    // Fallback to give normal behavior
    *first_slice = 0;
    *last_slice = ws->p;
    *zidx = 0;
    ws->skip = false;
#endif
}

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
