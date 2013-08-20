#include "SE_direct.h"
#include <math.h>

void SE3P_direct_real(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict x, 
		      const double* restrict q, int N,
		      const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
#ifdef _OPENMP
#pragma omp parallel for private(rvec, qn, p, r)
#endif
    for(int m=0; m<nidx; m++)
    {
	p = 0;
	for(int n=0; n<N; n++)
	{
	    rvec[0] = x[idx[m]    ]-x[n    ];
	    rvec[1] = x[idx[m]+N  ]-x[n+N  ];
	    rvec[2] = x[idx[m]+2*N]-x[n+2*N];
	    qn = q[n];

	    for(int p0 = -opt.layers; p0<=opt.layers; p0++)
		for(int p1 = -opt.layers; p1<=opt.layers; p1++)
		    for(int p2 = -opt.layers; p2<=opt.layers; p2++)
		    {
			if(idx[m] == n && p2 == 0 && p1 == 0 && p0 == 0)
			    continue;
			
			r = sqrt((rvec[0]+p0*opt.box[0])*
				 (rvec[0]+p0*opt.box[0])+
				 (rvec[1]+p1*opt.box[1])*
				 (rvec[1]+p1*opt.box[1])+
				 (rvec[2]+p2*opt.box[2])*
				 (rvec[2]+p2*opt.box[2]));
			
			p += qn*erfc(opt.xi*r)/r;
		    }
	}
	phi[m] += p;
    }
}

void SE3P_direct_real_rc(double* restrict phi, 
			 const int* restrict idx, int nidx,
			 const double* restrict x, 
			 const double* restrict q, int N,
			 const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
#ifdef _OPENMP
#pragma omp parallel for private(rvec, qn, p, r)
#endif
    for(int m=0; m<nidx; m++)
    {
	p = 0;
	for(int n=0; n<N; n++)
	{
	    rvec[0] = x[idx[m]    ]-x[n    ];
	    rvec[1] = x[idx[m]+N  ]-x[n+N  ];
	    rvec[2] = x[idx[m]+2*N]-x[n+2*N];
	    qn = q[n];

	    for(int p0 = -opt.layers; p0<=opt.layers; p0++)
		for(int p1 = -opt.layers; p1<=opt.layers; p1++)
		    for(int p2 = -opt.layers; p2<=opt.layers; p2++)
		    {
			if(idx[m] == n && p2 == 0 && p1 == 0 && p0 == 0)
			    continue;
			
			r = sqrt((rvec[0]+p0*opt.box[0])*
				 (rvec[0]+p0*opt.box[0])+
				 (rvec[1]+p1*opt.box[1])*
				 (rvec[1]+p1*opt.box[1])+
				 (rvec[2]+p2*opt.box[2])*
				 (rvec[2]+p2*opt.box[2]));

			if(r > opt.rc) continue;

			p += qn*erfc(opt.xi*r)/r;
		    }
	}
	phi[m] += p;
    }
}

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
	for(int j0 = -opt.layers; j0<=opt.layers; j0++)
	    for(int j1 = -opt.layers; j1<=opt.layers; j1++)
		for(int j2 = -opt.layers; j2<=opt.layers; j2++)
		{
		    if(j0 == 0 && j1 == 0 && j2==0)
			continue;

		    k[0] = 2*PI*j0/opt.box[0];
		    k[1] = 2*PI*j1/opt.box[1];
		    k[2] = 2*PI*j2/opt.box[2];
		    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
	
		    z=0;
		    for(int n = 0; n<N; n++)
			z += q[n]*cos(-(k[0]*(x[idx[m]]-x[n])+
					k[1]*(x[idx[m]+N  ]-x[n+N])+
					k[2]*(x[idx[m]+2*N]-x[n+2*N])));
	
		    p += z*exp(-k2/(4*opt.xi*opt.xi))/k2;
		}
	phi[m] += c*p;
    }
}

void SE3P_direct_self(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict q, int N, 
		      const ewald_opts opt)
{
    double c = 2*opt.xi/sqrt(PI);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m=0; m<nidx; m++)
	phi[m] -= c*q[idx[m]];
}
