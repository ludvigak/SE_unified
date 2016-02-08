#include "math.h"
#include "rotlet_direct.h"

#ifdef _OPENMP
#include "omp.h"
#endif

void rotlet_direct_rsrc(double* restrict u, 
			const double* restrict xt, 
			const int Nt,
			const double* restrict x, 
			const double* restrict f, 
			const int N,
			const ewald_opts opt)
{
    double r[3];
    double xm[3];
    int i1, i2, i3, m, n;
    const int nbox = 1;
    double rc2 = opt.rc * opt.rc;
    double xi = opt.xi;

#ifdef _OPENMP
#pragma omp parallel for \
    private(r,xm,i1,i2,i3,m,n) \
    default(shared)
#endif
    for(m=0; m<Nt; m++) // for all evaluation points
    {
	double um[3] = {0.0, 0.0, 0.0};

	xm[0] = xt[m     ];
	xm[1] = xt[m+Nt  ];
	xm[2] = xt[m+2*Nt];

	for(n=0; n<N; n++) // for all particles
	{
	    double xmn[3] = {xm[0]-x[n    ],
			     xm[1]-x[n+  N],
			     xm[2]-x[n+2*N]};
	    double f0 = f[n];
	    double f1 = f[n+N];
	    double f2 = f[n+2*N];
	    for(i1 = -nbox; i1<=nbox; i1++) // image boxes
		for(i2 = -nbox; i2<=nbox; i2++)
		    for(i3 = -nbox; i3<=nbox; i3++)
		    {
			// Assuming that r != 0 in home box
			r[0] = xmn[0]+opt.box[0]*i1;
			r[1] = xmn[1]+opt.box[1]*i2;
			r[2] = xmn[2]+opt.box[2]*i3;
			double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
			if(r2 > rc2)
			    continue; // skip outside rc
			double rnorm = sqrt(r2);
			double rxi = rnorm*xi;
			double A = (erfc(rxi)/rnorm + 
				    2*xi*exp(-rxi*rxi)/sqrt(PI) 
				    ) / r2;            
			um[0] += A*(f1*r[2] - f2*r[1]);
			um[1] += A*(f2*r[0] - f0*r[2]);
			um[2] += A*(f0*r[1] - f1*r[0]);
		    }
	}
	u[m     ] = um[0];
	u[m+Nt  ] = um[1];
	u[m+2*Nt] = um[2];
    }
}


