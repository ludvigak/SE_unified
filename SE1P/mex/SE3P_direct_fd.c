#include "mex.h"
#include "SE_direct.h"

#define IDX prhs[0]
#define X   prhs[1] // Source locations
#define Q   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define PHI plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

/* common option-unpacking */
void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    if(opt->xi==0)
      mexErrMsgTxt("xi cannot be zero");
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
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for(int i=0; i<N; i++)
	idx[i] = (int)idx_d[i] - 1;
}

double
heaviside(int x)
{
   if(x==0)
	return .5;
   else
	return 1.;
}

double 
hs(int x, int y, int z)
{
   double h=heaviside(x)*heaviside(y)*heaviside(z);
   if(h==1)
	return 8.;
   else if(h==.5)
	return 4.;
   else
	return 2;
}


#ifdef FORCE
void SE3P_direct_fd(double* restrict force,
		    const int* restrict idx, int nidx,
		    const double* restrict x,
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
    double c = 2.*PI/(opt.box[0]*opt.box[1]*opt.box[2]);
    double xi2 = opt.xi*opt.xi;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m=0; m<nidx; m++)
    {
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
	double f[3] = {0,0,0};
        for(int j0 = -opt.layers; j0<=opt.layers; j0++)
            for(int j1 = -opt.layers; j1<=opt.layers; j1++)
                for(int j2 = -opt.layers; j2<=opt.layers; j2++)
                {
                    if(j0 == 0 && j1 == 0 && j2==0)
                        continue;
		    double k[3];
                    k[0] = 2*PI*j0/opt.box[0];
                    k[1] = 2*PI*j1/opt.box[1];
                    k[2] = 2*PI*j2/opt.box[2];
                    double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                    double z=0;
                    for(int n = 0; n<N; n++){
                        z += q[n]*sin(k[0]*(x[n    ]-xm[0])+
                                      k[1]*(x[n+  N]-xm[1])+
                                      k[2]*(x[n+2*N]-xm[2]));
                     }

                    double h = 1.;

                    f[0] += h*z*exp(-k2/(4.*xi2))/k2*k[0];
                    f[1] += h*z*exp(-k2/(4.*xi2))/k2*k[1];
                    f[2] += h*z*exp(-k2/(4.*xi2))/k2*k[2];
                }

        force[idx[m]    ] = c*f[0];
        force[idx[m]+  N] = c*f[1];
        force[idx[m]+2*N] = c*f[2];
    }
}
#else
void SE3P_direct_fd(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
    double k[3];
    double k2, z, p;
    double c = 4*PI/(opt.box[0]*opt.box[1]*opt.box[2]);

#ifdef _OPENMP
#pragma omp parallel for private(k,k2,z,p)
#endif
    for(int m=0; m<nidx; m++)
    {
	p = 0;
	for(int j0 = 0; j0<=opt.layers; j0++)
	    for(int j1 = 0; j1<=opt.layers; j1++)
		for(int j2 = 0; j2<=opt.layers; j2++)
		{
		    if(j0 == 0 && j1 == 0 && j2==0)
			continue;
		    k[0] = 2*PI*j0/opt.box[0];
		    k[1] = 2*PI*j1/opt.box[1];
		    k[2] = 2*PI*j2/opt.box[2];
		    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
	
		    z=0;
		    for(int n = 0; n<N; n++)
                        z += q[n]*cos(k[0]*(x[idx[m]]-x[n]))*
                                     cos(k[1]*(x[idx[m]+N]-x[n+N]))*
                                     cos(k[2]*(x[idx[m]+2*N]-x[n+2*N]));

		    double h = hs(j0,j1,j2);
		    p += h*z*exp(-k2/(4*opt.xi*opt.xi))/k2;
		}
	phi[m] += c*p;
    }
}
#endif

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

#ifndef FORCE
    PHI = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
    double* restrict phi = mxGetPr(PHI);
#else 
    /* This is to allocate 3 vectors for the force. 
     * (FIXME) Note that the variable is still called PHI.*/
    PHI = mxCreateDoubleMatrix(num_eval, 3, mxREAL);
    double* restrict phi = mxGetPr(PHI);
#endif

    ewald_opts opt;
    unpack_opt(&opt, OPT);

    if(VERBOSE)
    {
	mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","FD1P",N,num_eval);
	mexPrintf("xi = %.2f, layers=%d\n", opt.xi,opt.layers);
    }

    // call kernel
    SE3P_direct_fd(phi, idx, num_eval, x, q, N, opt);
    mxFree(idx);
}
