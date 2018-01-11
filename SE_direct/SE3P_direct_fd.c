#include "math.h"
#include "mex.h"

#define IDX prhs[0]
#define X   prhs[1] // Source locations
#define Q   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define PHI plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define PI 3.141592653589793

typedef struct 
{
    double box[3];
    double xi;
    int layers; 
    double rc; 
} ewald_opts;

void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));

    opt->box[0] = box[0];
    opt->box[1] = box[1];
    opt->box[2] = box[2];

    // layers: mandatory for ewald sums that are truncated 
    const mxArray* mx_layers = mxGetField(mx_opt,0,"layers");
    if(mx_layers)
	opt->layers = (int)mxGetScalar(mx_layers);
    else
	opt->layers = -1;

    // rc: mandatory for short-range real sum 
    const mxArray* mx_rc = mxGetField(mx_opt,0,"real_cutoff");
    if(mx_rc)
	opt->rc = mxGetScalar(mx_rc);
    else
	opt->rc = -1;
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
}


/* no input checking is done */
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    // input dims
    const int N = mxGetM(X);
    
    const int num_eval = mxGetN(IDX); // FIXME: indices assumed to be row vec
    const double* idx_d = mxGetPr(IDX);
    int* idx = mxMalloc(num_eval*sizeof(int));
    index_translation(idx, idx_d, num_eval);
 
    const double* x = mxGetPr(X);
    const double* q = mxGetPr(Q);

    PHI = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
    double* restrict phi = mxGetPr(PHI);

    ewald_opts opt;
    unpack_opt(&opt, OPT);

    if(VERBOSE)
    {
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","FS3P",N,num_eval);
	mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
		  opt.xi,opt.rc,opt.layers);
    }
    // call kernel
    double k[3];
    double k2, z, p;
    double c = 4*PI/(opt.box[0]*opt.box[1]*opt.box[2]);
    double fac[3] = {2.*PI/opt.box[0], 2.*PI/opt.box[1], 2.*PI/opt.box[2]};
    double a = 1.0/(4*opt.xi*opt.xi);
#ifdef _OPENMP
#pragma omp parallel for private(k,k2,z,p)
#endif
    for(int m=0; m<num_eval; m++)
    {
	p = 0;
	double xm[3] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
	for(int j0 = -opt.layers; j0<=opt.layers; j0++)
	    for(int j1 = -opt.layers; j1<=opt.layers; j1++)
		for(int j2 = -opt.layers; j2<=opt.layers; j2++)
		{
		    if(j0 == 0 && j1 == 0 && j2==0)
			continue;

		    k[0] = fac[0]*j0;
		    k[1] = fac[1]*j1;
		    k[2] = fac[2]*j2;
		    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
	
		    z=0;
		    for(int n = 0; n<N; n++) {
		      double tmp = -(k[0]*(xm[0]-x[n]    )+
				     k[1]*(xm[1]-x[n+N]  )+
				     k[2]*(xm[2]-x[n+2*N])
				     );
		      z += q[n]*cos(tmp);
		    }
	
		    p += z*exp(-a*k2)/k2;
		}
	phi[m] += c*p;
    }    
    mxFree(idx);
}
