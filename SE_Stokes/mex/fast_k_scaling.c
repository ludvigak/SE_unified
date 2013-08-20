#include "mex.h"
#include "math.h"

#define G1   prhs[0]
#define G2   prhs[1]
#define G3   prhs[2]
#define XI   prhs[3] // Ewald parameters
#define L    prhs[4] // Box size
#define ETA  prhs[5]

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
    dims = mxGetDimensions(G1);
    //const int N = dims[0];
    //const int kmax = (N-1)/2; // NO CHECKING IS DONE THAT N IS ODD!!!
    
    const int N0 = dims[0];
	const int N1 = dims[1];
	const int N2 = dims[2];
    const int k0max = kmax(N0); // max wave number & index of zero wave
	const int k1max = kmax(N1);
	const int k2max = kmax(N2);
    
    const double xi = mxGetScalar(XI);
    const double eta = mxGetScalar(ETA);
    const double c = (1-eta)/(4*xi*xi);
    const double* box = mxGetPr(L);

    // real part of input arrays
    const double* restrict g1r = mxGetPr(G1);
    const double* restrict g2r = mxGetPr(G2);
    const double* restrict g3r = mxGetPr(G3);

    // imaginary parts
    const double* restrict g1i = mxGetPi(G1);
    const double* restrict g2i = mxGetPi(G2);
    const double* restrict g3i = mxGetPi(G3);

    // output arrays, complex
    H1 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);
    H2 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);
    H3 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);

    double* h1r = mxGetPr(H1);
    double* h2r = mxGetPr(H2);
    double* h3r = mxGetPr(H3);

    double* h1i = mxGetPi(H1);
    double* h2i = mxGetPi(H2);
    double* h3i = mxGetPi(H3);

    // scrach vars
    double B[3][3];
    double G[3];
    double k[3];

    int i0, i1, i2, idx;
    double k2;
    double q;


    for(i2 = 0; i2<N2; i2++)                  
    {           
	k[2] = 2.0*PI*(i2-k2max)/box[2];
	for(i1 = 0; i1<N1; i1++)                  
	{
	    k[1] = 2.0*PI*(i1-k1max)/box[1];
	    for(i0 = 0; i0<N0; i0++)
	    {
		idx = i0 + i1*N0 + i2*N0*N1;
		k[0] = 2.0*PI*(i0-k0max)/box[0];
		
		if(i0 != k0max || i1 != k1max || i2 != k2max)       // exclude k=0
		{

		    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
		    q = exp(-c*k2);
		    op_BB(B,k,xi);

		    // real part
		    G[0] = g1r[idx];
		    G[1] = g2r[idx];
		    G[2] = g3r[idx];

		    h1r[idx] = q*( B[0][0]*G[0]+B[0][1]*G[1]+B[0][2]*G[2] );
		    h2r[idx] = q*( B[1][0]*G[0]+B[1][1]*G[1]+B[1][2]*G[2] );
		    h3r[idx] = q*( B[2][0]*G[0]+B[2][1]*G[1]+B[2][2]*G[2] );

		    // imaginary part
		    G[0] = g1i[idx];
		    G[1] = g2i[idx];
		    G[2] = g3i[idx];

		    h1i[idx] = q*( B[0][0]*G[0]+B[0][1]*G[1]+B[0][2]*G[2] );
		    h2i[idx] = q*( B[1][0]*G[0]+B[1][1]*G[1]+B[1][2]*G[2] );
		    h3i[idx] = q*( B[2][0]*G[0]+B[2][1]*G[1]+B[2][2]*G[2] );

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
