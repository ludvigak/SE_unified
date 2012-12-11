#include <math.h>
#include "SE_Stokes_direct.h"

#ifdef HASIMOTO
#include "hasimoto_decomp.h"
#elif BEENAKKER
#include "beenakker_decomp.h"
#else
#error "Must provide -D<decomposition> to compiler"
#endif

void SE3P_Stokes_direct_real(double* restrict u, 
			     const int* restrict idx, int nidx,
			     const double* restrict x, 
			     const double* restrict f, int N,
			     const ewald_opts opt)
{
    const int nbox = opt.layers;

    double r[3];
    double xm[3];
    double A[3][3];
    int i1, i2, i3, m, n;

    for(m=0; m<nidx; m++)                          // for all evaluation points
    {
	u[m       ] = 0;
	u[m+nidx  ] = 0;
	u[m+2*nidx] = 0;

	xm[0] = x[idx[m]    ];             // indirect indexing OK in outer loop
	xm[1] = x[idx[m]+N  ];
	xm[2] = x[idx[m]+2*N];

	for(i1 = -nbox; i1<=nbox; i1++)                           // image boxes
	    for(i2 = -nbox; i2<=nbox; i2++)
		for(i3 = -nbox; i3<=nbox; i3++)
		{
		    for(n=0; n<N; n++)                      // for all particles
		    {
			if(i1==0 && i2==0 && i3==0 && n==idx[m])    // skip self
			    continue;

			r[0] = xm[0]-x[n    ]+opt.box[0]*i1;
			r[1] = xm[1]-x[n+  N]+opt.box[1]*i2;
			r[2] = xm[2]-x[n+2*N]+opt.box[2]*i3;

			op_A(A,r,opt.xi);                            // u += A*f
			u[m       ] += 
			    A[0][0]*f[n]+A[0][1]*f[n+N]+A[0][2]*f[n+2*N];
			u[m+nidx  ] += 
			    A[1][0]*f[n]+A[1][1]*f[n+N]+A[1][2]*f[n+2*N];
			u[m+2*nidx] +=
			    A[2][0]*f[n]+A[2][1]*f[n+N]+A[2][2]*f[n+2*N];
		    }
		}
    }
}

void SE3P_Stokes_direct_real_rc(double* restrict u, 
				const int* restrict idx, int nidx,
				const double* restrict x, 
				const double* restrict f, int N,
				const ewald_opts opt)
{
    const int nbox = opt.layers;

    double r[3];
    double xm[3];
    double A[3][3];
    int i1, i2, i3, m, n;

    for(m=0; m<nidx; m++)                          // for all evaluation points
    {
	u[m       ] = 0;
	u[m+nidx  ] = 0;
	u[m+2*nidx] = 0;

	xm[0] = x[idx[m]    ];             // indirect indexing OK in outer loop
	xm[1] = x[idx[m]+N  ];
	xm[2] = x[idx[m]+2*N];

	for(i1 = -nbox; i1<=nbox; i1++)                           // image boxes
	    for(i2 = -nbox; i2<=nbox; i2++)
		for(i3 = -nbox; i3<=nbox; i3++)
		{
		    for(n=0; n<N; n++)                      // for all particles
		    {
			if(i1==0 && i2==0 && i3==0 && n==idx[m])    // skip self
			    continue;

			r[0] = xm[0]-x[n    ]+opt.box[0]*i1;
			r[1] = xm[1]-x[n+  N]+opt.box[1]*i2;
			r[2] = xm[2]-x[n+2*N]+opt.box[2]*i3;

			if(sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) > opt.rc)
			    continue;                         // skip outside rc

			op_A(A,r,opt.xi);                            // u += A*f
			u[m       ] += 
			    A[0][0]*f[n]+A[0][1]*f[n+N]+A[0][2]*f[n+2*N];
			u[m+nidx  ] += 
			    A[1][0]*f[n]+A[1][1]*f[n+N]+A[1][2]*f[n+2*N];
			u[m+2*nidx] +=
			    A[2][0]*f[n]+A[2][1]*f[n+N]+A[2][2]*f[n+2*N];
		    }
		}
    }
}

void SE3P_Stokes_direct_fd(double* restrict u, 
			   const int* restrict idx, int nidx,
			   const double* restrict x, 
			   const double* restrict f, int N,
			   const ewald_opts opt)
{
    double B[3][3];
    double z[3];
    double k[3];
    double xm[3];
    int i1, i2, i3, m, n;
    double q, k_dot_r;

    const double vol = opt.box[0]*opt.box[1]*opt.box[2];
    const int kmax=opt.layers;
    const double kc0 = 2.0*PI/opt.box[0];
    const double kc1 = 2.0*PI/opt.box[1];
    const double kc2 = 2.0*PI/opt.box[2];

    for(m=0; m<nidx; m++)                          // for all evaluation points
    {
	u[m       ] = 0;
	u[m+  nidx] = 0;
	u[m+2*nidx] = 0;

	xm[0] = x[idx[m]    ];             // indirect indexing OK in outer loop
	xm[1] = x[idx[m]+N  ];
	xm[2] = x[idx[m]+2*N];

	for(i1 = -kmax; i1<=kmax; i1++)                      // for k-space cube
	    for(i2 = -kmax; i2<=kmax; i2++)
		for(i3 = -kmax; i3<=kmax; i3++)
		{
		    if(i3 != 0 || i2 != 0 || i1 != 0)             // exclude k=0
		    {
			z[0] = 0; z[1] = 0; z[2] = 0;
			k[0] = kc0*i1;
			k[1] = kc1*i2; 
			k[2] = kc2*i3;

			for(n=0; n<N; n++)                  // for all particles
			{
			    k_dot_r = 
				k[0]*(xm[0]-x[n    ])+ 
				k[1]*(xm[1]-x[n+N  ])+ 
				k[2]*(xm[2]-x[n+2*N]);
			    q = cos(k_dot_r);
			    z[0] += q*f[n    ];
			    z[1] += q*f[n+  N];
			    z[2] += q*f[n+2*N];
			}
			op_B(B,k,opt.xi);                      // multiplication
			u[m       ] += B[0][0]*z[0]+B[0][1]*z[1]+B[0][2]*z[2];
			u[m+nidx  ] += B[1][0]*z[0]+B[1][1]*z[1]+B[1][2]*z[2];
			u[m+2*nidx] += B[2][0]*z[0]+B[2][1]*z[1]+B[2][2]*z[2];
		    }
		}
	u[m       ] /= vol;
	u[m+  nidx] /= vol;
	u[m+2*nidx] /= vol;
    }
}

void SE3P_Stokes_direct_self(double* restrict u, 
			     const int* restrict idx, int nidx,
			     const double* restrict f, int N, 
			     const ewald_opts opt)
{
    double c = self_coeff(opt.xi);
    for(int m=0; m<nidx; m++)
    {
	u[m       ] = c*f[idx[m]    ];
	u[m+  nidx] = c*f[idx[m]+N  ];
	u[m+2*nidx] = c*f[idx[m]+2*N];
    }
}
