#include "mex.h"
#include "math.h"

#define G1   prhs[0]
#define G2   prhs[1]
#define G3   prhs[2]
#define GR   prhs[3]
#define XI   prhs[4] // Ewald parameters
#define L    prhs[5] // Box size
#define ETA  prhs[6]

#define H1   plhs[0]  // Output
#define H2   plhs[1]  // Output
#define H3   plhs[2]  // Output

#define PI 3.141592653589793
#define NDIMS 3

// Hasimoto decomposition for stokeslet
double inline op_poly(double k2, double xi)
{
    return 1 + k2/(4*xi*xi);
}
void inline op_tensor(double H[3], double G[3], double k[3], double k2)
{
    double kdotg = k[0]*G[0] + k[1]*G[1] + k[2]*G[2];
    H[0] = (k2*G[0] - kdotg*k[0]);
    H[1] = (k2*G[1] - kdotg*k[1]);
    H[2] = (k2*G[2] - kdotg*k[2]);
}


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

    const mwSize *dims;
    dims = mxGetDimensions(G1);
    
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
    const double* restrict gR = mxGetPr(GR);

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

    double* restrict h1r = mxGetPr(H1);
    double* restrict h2r = mxGetPr(H2);
    double* restrict h3r = mxGetPr(H3);

    double* restrict h1i = mxGetPi(H1);
    double* restrict h2i = mxGetPi(H2);
    double* restrict h3i = mxGetPi(H3);

    // scratch vars
    double G[3];
    double H[3];
    double k[3];
#pragma omp parallel for private(G,H,k)
    for(int i2 = 0; i2<N2; i2++)                  
    {           
	k[2] = 2.0*PI*(i2-k2max)/box[2];
	for(int i1 = 0; i1<N1; i1++)                  
	{
	    k[1] = 2.0*PI*(i1-k1max)/box[1];
	    for(int i0 = 0; i0<N0; i0++)
	    {
		k[0] = 2.0*PI*(i0-k0max)/box[0];		
		int idx = i0 + i1*N0 + i2*N0*N1;
		double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
		double q = exp(-c*k2);
		double B = q * gR[idx] * op_poly(k2, xi);
		// real part
		G[0] = g1r[idx];
		G[1] = g2r[idx];
		G[2] = g3r[idx];
		op_tensor(H, G, k, k2);
		h1r[idx] = B*H[0];
		h2r[idx] = B*H[1];
		h3r[idx] = B*H[2];
		// imaginary part
		G[0] = g1i[idx];
		G[1] = g2i[idx];
		G[2] = g3i[idx];
		op_tensor(H, G, k, k2);
		h1i[idx] = B*H[0];
		h2i[idx] = B*H[1];
		h3i[idx] = B*H[2];
	    }
	}
    }
}
