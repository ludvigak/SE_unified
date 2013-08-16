#include "stresslet_real_rc.h"

static double randnum(double, double);

int main() {
     int N = 100000;
     double rc = 0.1;
     double xi = 5.0;
     double L = 1;
    // Setup system
    double* restrict x = malloc(3*N*sizeof(double));
    double* restrict nvec = malloc(3*N*sizeof(double));
    double* restrict fvec = malloc(3*N*sizeof(double));

    for(int i=0;i<3*N;i++)
    {
	// Cluster particels
	x[i] = L*randnum(0,1);
	nvec[i] = randnum(0,1);
	fvec[i] = randnum(0,1);
    }

    double* restrict box = malloc(3*sizeof(int));
    box[0] = L;
    box[1] = L;
    box[2] = L;

    int numel;
    int* restrict row;
    int* restrict col;
    int* restrict idx_in_array;
    int* restrict buck_size;
    double* restrict val[3][3];

    // Get the triplets of points in the matrices
    /* get_rs_triplets ( x, nvec, N,  */
    /* 		       box, xi, rc, 1, */
    /* 		       &row, &col, val, &buck_size, &idx_in_array, &numel */
    /* 		       ); */

    // Compute matrix free
    double* restrict phi; 
    compute_rsrc_direct (x, nvec, fvec, N, box, xi, rc, &phi);

    double* restrict phi_new; 
    compute_rsrc_direct_new (x, nvec, fvec, N, box, xi, rc, &phi_new);

    double maxerr = 0.0;
    double maxval = 0.0;
    for (int i=0; i<N; i++)
    {
	if (abs(phi[i])>maxval)
	    maxval = abs(phi[i]);
	double err = abs(phi[i]-phi_new[i]);
	if(err > maxerr)
	    maxerr = err;
    }
    printf("maxerr_rel = %g\n maxval=%g\n", maxerr, maxval);

    return 0;
}

static double randnum(double min, double L)
{
    double q = ( (double) rand() )/RAND_MAX;
    return L*q+min;
}
