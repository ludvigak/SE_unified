#include "mex.h"
#include "math.h"

#define X    prhs[0] // Source locations
#define Q    prhs[1] // Source strengths
#define N    prhs[2] // First N elements to compute at
#define NBH  prhs[3] // Neighbor list (column indices in x and f)
#define LEN  prhs[4]
#define XI   prhs[5] // Ewald parameters

#define PHI plhs[0]  // Output

#define PI 3.141592653589793

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    if(VERBOSE)
	mexPrintf("Ewald kernel (real, nbh)\n");
    
    const double* restrict x   = mxGetPr(X);
    const double* restrict q   = mxGetPr(Q);
    const int n = (int) mxGetScalar(N);
    const int* restrict nbh = (int*) mxGetPr(NBH); /* DANGER! */
    const int* restrict len = (int*) mxGetPr(LEN);
    const double xi = mxGetScalar(XI);

    PHI = mxCreateDoubleMatrix(n, 1, mxREAL);
    double* restrict phi = mxGetPr(PHI);

    double r[3], d;

    int i,j,idx,k=0,i0=0;
    for(i=0; i<n; i++) /* for each x */
    {
	/*
	  THIS LOOP IS NOT HPC!
	  i)  Inner-loop brach makes SSE-vectorization impossible. It also
	  causes a stall (probably).
	  ii) Indirect indexing x[nbh[k]] will ruin hw prefetch and trash 
	  the cache. 
	  
	  UPDATE:
	  Runtime is dominated by the cost of computing 'erfc',
	  masking the cache-miss latency issue.
	*/
	phi[i]=0;
	for(j=0; j<len[i0]; j++) /* for each neighbor */
	{
	    idx = nbh[k++]; /* neighbor index */

	    if(i==idx)
	    	continue;

	    r[0] = x[i    ] - x[idx    ];
	    r[1] = x[i+  n] - x[idx+  n];
	    r[2] = x[i+2*n] - x[idx+2*n];

	    d = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	    
	    phi[i] += q[idx]*erfc(xi*d)/d;
	}
	i0++;
    }
}
