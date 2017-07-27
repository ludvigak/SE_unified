#include "mex.h"
#include "math.h"

#define G33  prhs[0]
#define GR   prhs[1]
#define XI   prhs[2] // Ewald parameters
#define L    prhs[3] // Box size
#define ETA  prhs[4]

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
void inline op_tensor(double H[3], double G[3][3], double k[3], double k2, double B)
{
    for (int j=0; j<3; j++)	 
    {
	double tmp1 = 0;
	double tmp2 = 0;
	for (int l=0; l<3; l++)
	{
	    tmp1 = tmp1 + k[j]*G[l][l] + k[l]*(G[j][l] + G[l][j]);
	    tmp2 = tmp2 - k[l]*( k[0]*G[l][0] + k[1]*G[l][1] + k[2]*G[l][2]);
	}
	H[j] = B*(k2*tmp1 + 2*k[j]*tmp2);
    }
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

void inline copy_in(double G[3][3], const double* restrict g[3][3], size_t idx)
{
    for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)	 
	    G[i][j] = g[i][j][idx];
}

void inline copy_out(double H[3], double* restrict h[3], size_t idx)
{
    for (int i=0; i<3; i++)
	h[i][idx] = H[i];
}

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{

    // Input G33 is 3x3 cell array, get input arrays
    const double* restrict gr[3][3];
    const double* restrict gi[3][3];

    for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)	    
	{
	    mxArray* cell = mxGetCell(G33, i*3+j);
	    if(mxIsComplex(cell)==false)
		mexErrMsgTxt("Inputs must be complex!");
	    gr[i][j] = mxGetPr(cell);
	    gi[i][j] = mxGetPi(cell);	    
	}
    
    // Setup sizes
    const mwSize *dims;
    dims = mxGetDimensions(GR);
    
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

    // output arrays, complex
    H1 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);
    H2 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);
    H3 = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS,mxCOMPLEX);
    double* restrict hr[3];
    double* restrict hi[3];    
    hr[0] = mxGetPr(H1);
    hr[1] = mxGetPr(H2);
    hr[2] = mxGetPr(H3);
    hi[0] = mxGetPi(H1);
    hi[1] = mxGetPi(H2);
    hi[2] = mxGetPi(H3);

    // scratch vars
    double G[3][3];
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
		copy_in(G, gi, idx);
		op_tensor(H, G, k, k2, B);
		copy_out(H, hr, idx);
		// imaginary part
		copy_in(G, gr, idx);
		op_tensor(H, G, k, k2, -B);
		copy_out(H, hi, idx);
	    }
	}
    }
}
