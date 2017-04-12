#include "mex.h"
#include "mex_compat.h"
#include "math.h"
#include "cell_list.h"

#ifdef INTEL_MKL
#include "mkl.h"
#endif


#define X   prhs[0] // Source locations
#define F   prhs[1] // Source strengths
#define RC  prhs[2] // cutoff
#define XI  prhs[3] // Ewald Param
#define P   prhs[4] // Periodic wrap
#define BOX prhs[5] // domain size

#define U   plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

#ifdef _OPENMP
#define CRITICAL _Pragma("omp critical")
#else
#define CRITICAL
#endif

#define PI 3.141592653589793

inline double norm2(double * a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

// Compute buffer stuff
#define BUF_SIZE 256
typedef struct {
    int n;
    int idx_t[BUF_SIZE];
    double rvec[3*BUF_SIZE];
    double rsq[BUF_SIZE];
} ComputeBuffer;

static void buffer_push(ComputeBuffer* buffer, int idx_t, double* rvec, double rsq)
{
    int n = buffer->n;
    buffer->idx_t[n] = idx_t;
    for(int i=0; i<3; i++)
	buffer->rvec[3*n + i] = rvec[i];
    buffer->rsq[n] = rsq;
    buffer->n = n + 1;
}

static void empty_buffer(ComputeBuffer* buffer, 
			 const double* restrict x, 
			 const double* restrict f, 
			 double* restrict u, 
			 int idx_s,
			 double xi)
{
    int N = buffer->n;
    int idx_t;
    double fs, ft;
    double us =  0;
    fs = f[idx_s];
    
    // Do what we can to help the compiler vectorize exp and erfc, if possible
    const double* restrict r2 = buffer->rsq;
    double c1[BUF_SIZE];
#ifdef INTEL_MKL
    double r[BUF_SIZE];
    double xir[BUF_SIZE];
    double xi2r2[BUF_SIZE];
    for (int n=0; n<N; n++)
    {
	r[n] = sqrt(r2[n]);
	xir[n] = xi*r[n];
	xi2r2[n] = -xi2*r2[n];	
    }
    double erfc_vec[BUF_SIZE];
    double exp_vec[BUF_SIZE];
    vdErfc(N, xir, erfc_vec);
    vdExp(N, xi2r2, exp_vec);
    for (int n=0; n<N; n++)
    {
	double  xiexp = xi*exp_vec[n];
	c1[n] = erfc_vec[n] / r[n];
    }    
#else
    for (int n=0; n<N; n++)
    {     
      double r = sqrt(r2[n]);
      c1[n] = erfc(xi*r) / r;
    }
#endif
    // Compute interactions
    for (int n=0; n<N; n++)
    {
	idx_t = buffer->idx_t[n];
	ft = f[idx_t];
	u[idx_t] += fs*c1[n];
	us += ft*c1[n];
    }
    u[idx_s] += us;
    buffer->n = 0;
}

// Entry point

void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{

    // input target
    const int N = mxGetN(X);
    const double xi = (double) mxGetScalar(XI);
    const double rc = (double) mxGetScalar(RC);
    const double rcsq = rc*rc;
    const double p = (double) mxGetScalar(P);
    const double* x = mxGetPr(X);
    const double* f = mxGetPr(F);
    const double* box = mxGetPr(BOX);
  
    // output
    U = mxCreateDoubleMatrix(N, 1, mxREAL);
    double* restrict u_out = mxGetPr(U);

    // Setup cell list variables
    int ncell[3];
    int* restrict cell_list;
    int* restrict cell_idx;
    double rn;
    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // Build cell list
    build_cell_list(x, N, box, rc, &rn, ncell, &cell_list, &cell_idx);
  

#ifdef _OPENMP
#pragma omp parallel
#endif
    { // Begin parallel section    
	// Create a thread-local compute buffer
	ComputeBuffer buffer;
	// Setup local output
	double* restrict u;
	CRITICAL {
	  u = __MALLOC(N*sizeof(double));
	}
	for(int i=0;i<N;i++)
	    u[i] = 0.0;

	// Main loop
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
	for (int idx_s=0; idx_s<N; idx_s++)  
	{
	    double xs[3];
	    int home_cell[3], icell[3];

	    for(int i=0; i<3; i++)
	    {
		// Source point
		xs[i] = x[idx_s*3 + i];
		// Determine home cell
		home_cell[i] = xs[i]/rn;	
	    }

	    // Iterate through near cells (including home cell)
	    buffer.n = 0;
	    for(int ip=0; ip<27; ip++) 
	    {
		// Get neigh cell
		icell[0] = home_cell[0] + px[ip];
		icell[1] = home_cell[1] + py[ip];
		icell[2] = home_cell[2] + pz[ip];
				
		// Stop at boundaries
		int inside = 1;
		for(int j=0; j<3; j++)
		{
		  if (icell[j] < 0 || icell[j] == ncell[j])
			inside = 0;
		}
		if (!inside)
		    continue;
		int icell_idx = 
		    icell[0] +
		    icell[1]*ncell[0] + 
		    icell[2]*ncell[1]*ncell[0];	    
		// Go through cell list
		int cell_a = cell_idx[icell_idx];
		int cell_b = cell_idx[icell_idx+1];
		for(int point_idx=cell_a; point_idx<cell_b; point_idx++)
		{
		    int idx_t = cell_list[point_idx];
		    if (idx_s >= idx_t) 
			continue;
		    double rvec[3];
		    // periodic wrap
		    for (int j=-p; j<=p; j++)
		    {
		      double pshift[] = {0,0,j*box[2]};
		      for(int i=0; i<3; i++)
			rvec[i] = xs[i] - x[idx_t*3 + i] - pshift[i];
		      double r2 = norm2(rvec);
		      if (r2 > rcsq) 
			continue; 
		      buffer_push(&buffer, idx_t, rvec, r2);
		      if (buffer.n == BUF_SIZE)
			empty_buffer(&buffer, x, f, u, idx_s, xi);
		    }
		}
	    }
	    empty_buffer(&buffer, x, f, u, idx_s, xi);
	}
	// End of particle loop, collect results
	CRITICAL {
	    for(int i=0; i<N; i++)
		u_out[i] += u[i];
	}
	// free/malloc not thread safe under MEX
	CRITICAL {
	    __FREE(u);
	}
    }
}
