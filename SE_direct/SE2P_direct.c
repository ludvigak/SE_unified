#include "SE_direct.h"
#include <math.h>

#define __EXP_ARG_MAX 600

void SE2P_direct_real(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict x, 
		      const double* restrict q, int N,
		      const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
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
		{
		    if(idx[m] == n && p1 == 0 && p0 == 0)
			continue;
			
		    r = sqrt((rvec[0]+p0*opt.box[0])*
			     (rvec[0]+p0*opt.box[0])+
			     (rvec[1]+p1*opt.box[1])*
			     (rvec[1]+p1*opt.box[1])+
			     rvec[2]*rvec[2]);
			
		    p += qn*erfc(opt.xi*r)/r;
		}
	}
	phi[m] = p;
    }
}

void SE2P_direct_real_rc(double* restrict phi, 
			 const int* restrict idx, int nidx,
			 const double* restrict x, 
			 const double* restrict q, int N,
			 const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
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
		{
		    if(idx[m] == n && p1 == 0 && p0 == 0)
			continue;
			
		    r = sqrt((rvec[0]+p0*opt.box[0])*
			     (rvec[0]+p0*opt.box[0])+
			     (rvec[1]+p1*opt.box[1])*
			     (rvec[1]+p1*opt.box[1])+
			     rvec[2]*rvec[2]);
			
		    if(r > opt.rc) continue;

		    p += qn*erfc(opt.xi*r)/r;
		}
	}
	phi[m] = p;
    }
}

static inline double theta_plus(double z, double k, double xi)
{
    /* idea for a more stable form [LK] */
    /* exp( k*z + log( erfc(k/(2.0*xi) + xi*z) ) ); */

    if(k*z <  __EXP_ARG_MAX)
	return exp( k*z)*erfc(k/(2.0*xi) + xi*z);
    else 
	return 0.0;
}

static inline double theta_minus(double z, double k, double xi)
{
    /* exp(-k*z + log( erfc(k/(2.0*xi) - xi*z) ) ); */

    if(-k*z <  __EXP_ARG_MAX)
	return exp(-k*z)*erfc(k/(2.0*xi) - xi*z);
    else 
	return 0.0;
}

void SE2P_direct_fd(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
    double k[2], xm[3]; 
    double kn, k_dot_r, z, phi_m;
    double cm, cp;
    const double xi = opt.xi;

    for(int m=0; m<nidx; m++)
    {
	xm[0] = x[idx[m]    ];
	xm[1] = x[idx[m]+N  ];
	xm[2] = x[idx[m]+2*N];
	phi_m = 0;
	for(int n = 0; n<N; n++){
	    for(int j0 = -opt.layers; j0<=opt.layers; j0++)
		for(int j1 = -opt.layers; j1<=opt.layers; j1++)
		{
		    if(j0 == 0 && j1 == 0)
			continue;

		    k[0] = 2*PI*j0/opt.box[0];
		    k[1] = 2*PI*j1/opt.box[1];
		    kn = sqrt(k[0]*k[0] + k[1]*k[1]);
		    k_dot_r = k[0]*(xm[0]-x[n]) + k[1]*(xm[1]-x[n+N]);
		    z = xm[2]-x[n+2*N];
		    cp = theta_plus(z,kn,xi);
		    cm = theta_minus(z,kn,xi);

		    phi_m += q[n]*cos(k_dot_r)*(cm+cp)/kn;
		}
	}
	phi[m] = PI*phi_m/(opt.box[0]*opt.box[1]);
    }
}

void SE2P_direct_k0(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
    double z, zm, phi_m;
    for(int m=0; m<nidx; m++)
    {
	phi_m=0;
	zm = x[idx[m]+2*N];
	for(int n=0; n<N; n++)
	{
	    z = zm - x[n+2*N];
	    phi_m += q[n]*(exp(-opt.xi*opt.xi*z*z)/opt.xi + 
			   sqrt(PI)*z*erf(opt.xi*z));
	}
	phi[m] = -2*phi_m*sqrt(PI)/(opt.box[0]*opt.box[1]);
    }
}

void SE2P_direct_self(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict q, int N, 
		      const ewald_opts opt)
{
    double c = 2*opt.xi/sqrt(PI);
    for(int m=0; m<nidx; m++)
	phi[m] = -c*q[idx[m]];
}
