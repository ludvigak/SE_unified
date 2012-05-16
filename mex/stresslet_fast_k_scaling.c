#include "mex.h"
#include "math.h"

#define G    prhs[0]
#define XI   prhs[1] // Ewald parameters
#define L    prhs[2] // Box size
#define ETA  prhs[3]

#define H1   plhs[0]  // Output
#define H2   plhs[1]  // Output
#define H3   plhs[2]  // Output

#define PI 3.141592653589793
#define NDIMS 3

// select k-space scaling tensor
#ifdef HASIMOTO
#include "hasimoto_op_fd.h"
#elif BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

static inline int is_odd(int p)
{
    return p&1;
}

static int kmax(int N)
{
    if(is_odd(N))
    {
        return (N-1)/2;
    }
    else
    {
       return N/2;
    }
}


void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    if(VERBOSE)
	mexPrintf("Stokes Ewald FFT k-scaling: %s\n", OP_TAG);

    const mwSize *dims;
    dims = mxGetDimensions(G); // 9 x N0 x N1 x N2
    const int N0 = dims[1];
	const int N1 = dims[2];
	const int N2 = dims[3];
    const int k0max = kmax(N0); // max wave number & index of zero wave
	const int k1max = kmax(N1);
	const int k2max = kmax(N2);
    const double xi = mxGetScalar(XI);
    const double eta = mxGetScalar(ETA);
    const double c = (1-eta)/(4*xi*xi);
    const double* box = mxGetPr(L);

    // real part of input array
    const double* restrict gr = mxGetPr(G);

    // imaginary part
    const double* restrict gi = mxGetPi(G);

    // output arrays, complex
    H1 = mxCreateNumericArray(NDIMS,&dims[1],mxDOUBLE_CLASS,mxCOMPLEX);
    H2 = mxCreateNumericArray(NDIMS,&dims[1],mxDOUBLE_CLASS,mxCOMPLEX);
    H3 = mxCreateNumericArray(NDIMS,&dims[1],mxDOUBLE_CLASS,mxCOMPLEX);

    double* h1r = mxGetPr(H1);
    double* h2r = mxGetPr(H2);
    double* h3r = mxGetPr(H3);

    double* h1i = mxGetPi(H1);
    double* h2i = mxGetPi(H2);
    double* h3i = mxGetPi(H3);

    // scratch vars
    double B[27];
    double T[9];
    double k[3];

    int i0, i1, i2, idx, j;
    double k2;
    double q;
	double tmp;


    for(i2 = 0; i2<N2; i2++)                  
    {           
	k[2] = -2.0*PI*(i2-k2max)/box[2];
	for(i1 = 0; i1<N1; i1++)                  
	{
	    k[1] = -2.0*PI*(i1-k1max)/box[1];
	    for(i0 = 0; i0<N0; i0++)
	    {
		idx = i0 + i1*N0 + i2*N0*N1;
		k[0] = -2.0*PI*(i0-k0max)/box[0];
		
		if(i0 != k0max || i1 != k1max || i2 != k2max)       // exclude k=0
		{

		    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
		    q = exp(-c*k2);
		    op_BB(B,k,xi);

		    // real part of h
			for(j=0; j<9; j++)
			{
			    T[j] = gi[9*idx+j];
			}
			
			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j]*T[j];
			h1r[idx] = -q*tmp;

			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j+9]*T[j];
			h2r[idx] = -q*tmp;

			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j+18]*T[j];
			h3r[idx] = -q*tmp;

		    // imaginary part of h
			for(j=0; j<9; j++)
			{
			    T[j] = gr[9*idx+j];
			}
			
			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j]*T[j];
			h1i[idx] = q*tmp;

			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j+9]*T[j];
			h2i[idx] = q*tmp;

			tmp=0;
			for(j=0; j<9; j++)
				tmp += B[j+18]*T[j];
			h3i[idx] = q*tmp;

		}
		else // don't assume initalized to zero
		{
		    h1r[idx] = 0; h2r[idx] = 0; h3r[idx] = 0;
		    h1i[idx] = 0; h2i[idx] = 0; h3i[idx] = 0;
		}
	    }
	}
    }
}
