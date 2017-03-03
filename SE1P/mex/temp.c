#include <stdlib.h>
#include "SE_direct.h"
#include <math.h>
#include "mathint.h"

#include <omp.h>

void SE1P_direct(double* restrict phi, 
		 const int* restrict idx, int nidx,
		 const double* restrict x, 
		 const double* restrict q, int N,
		 const ewald_opts opt)
{
  double p;
#ifdef _OPENMP
#pragma omp parallel for private(p)
#endif
  for(int m=0; m<nidx; m++)
    {
      p=0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p2 = -opt.layers; p2<=opt.layers; p2++)
	    {
	      if(idx[m] == n && p2 == 0)
		continue;

	      double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      
	      p += qn/r;
	    }
	}
      phi[idx[m]] = p;
    }
}

void SE1P_direct_real(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict x, 
		      const double* restrict q, int N,
		      const ewald_opts opt)
{
    double rvec[3];
    double qn;
    double p, r;
#ifdef _OPENMP
#pragma omp parallel for private(p)
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

	    for(int p2 = -opt.layers; p2<=opt.layers; p2++)
		{
		    if(idx[m] == n && p2 == 0)
			continue;
			
		    r = sqrt((rvec[2]+p2*opt.box[2])*
			     (rvec[2]+p2*opt.box[2])+
			     rvec[0]*rvec[0]+
			     rvec[1]*rvec[1]);
			
		    p += qn*erfc(opt.xi*r)/r;
		}
	}
	phi[m] = p;
    }
}

void SE1P_direct_real_rc(double* restrict phi, 
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

	    for(int p2 = -opt.layers; p2<=opt.layers; p2++)
		{
		    if(idx[m] == n && p2 == 0)
			continue;
			
		    r = sqrt((rvec[2]+p2*opt.box[2])*
			     (rvec[2]+p2*opt.box[2])+
			     rvec[0]*rvec[0]+
			     rvec[1]*rvec[1]);
			
		    if(r > opt.rc) continue;

		    p += qn*erfc(opt.xi*r)/r;
		}
	}
	phi[m] = p;
    }
}

void SE1P_direct_fd(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
  double xm[3]; 
  double k3, k3z, z, phi_m, rho2, a, b, K0,qn;

  const double xi   = opt.xi;
  double xi2        = xi*xi;
  double TwoPiOverL = 2.*PI/opt.box[2];
  int rep;

  for(int m=0; m<nidx; m++)
    {
      xm[0] = x[idx[m]    ];
      xm[1] = x[idx[m]+N  ];
      xm[2] = x[idx[m]+2*N];
      phi_m = 0;
      for(int n = 0; n<N; n++)
	{
	  z   = xm[2]-x[n+2*N];
	  rho2= ( (xm[0]-x[n  ])*(xm[0]-x[n  ])+
	          (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  b   = rho2*xi2;
	  qn  = q[n];
	  for(int j2 = -opt.layers; j2<=opt.layers; j2++)
	    {
	      if(j2 == 0)
		continue;
	      
	      k3  = TwoPiOverL*j2;
	      k3z = -k3*z;
	      
	      a   = k3*k3/(4.*xi2);
	    
	      K0  = IncompBesselK0_Simpson(1e-15 ,&rep, a, b,0);
	      phi_m += qn*cos(k3z)*K0;
	    }
	}
      phi[m] = phi_m/(opt.box[2]);
    }
}

void SE1P_direct_fd0(double* restrict phi, 
		    const int* restrict idx, int nidx,
		    const double* restrict x, 
		    const double* restrict q, int N,
		    const ewald_opts opt)
{
  double xm[3],M; 
  double k3, k3z, z, phi_m, rho2, a, b, K0,qn;

  const double xi   = opt.xi;
  double xi2        = xi*xi;
  double TwoPiOverL = 2.*PI/opt.box[2];
  int xiL = xi*opt.box[2];
  
  // prepare for integration
//  M = 1000;
//  double *X = (double*) malloc(M*sizeof(double));  // quadrature points
//  double *W = (double*) malloc(M*sizeof(double));  // quadrature weights
/*
  int M1 = M/4;
  glwt_fast(M1,0.,.25,X,W);
  glwt_fast(M1,.25,.5,&X[M1],&W[M1]);
  glwt_fast(M1,.5,.75,&X[2*M1],&W[2*M1]);
  glwt_fast(M1,.75,1.,&X[3*M1],&W[3*M1]);
*/

  M         = 1000;                                  // Initial guess
  int M1    = MIN(700,MAX(xiL,200));                 
  int M2    = MIN(M-M1,300);
  // if M1+M2 is not perfect divisior of 4 then make it
  int di    = (M1+M2)%4; di = di%4;
  M2        = M2 + di;                               // remainder added to M2
  double v  = 1./(xi*opt.box[2]);	             // splitting point
  M         = M1+M2; // update the guess
  double *X = (double*) _mm_malloc(M*sizeof(double),32); // quadrature points
  double *W = (double*) _mm_malloc(M*sizeof(double),32); // quadrature weights

  glwt_fast(M1,0.,v,X,W);
  glwt_fast(M2,v,1.,&X[M1],&W[M1]);

  for(int m=0; m<nidx; m++)
    {
      xm[0] = x[idx[m]    ];
      xm[1] = x[idx[m]+N  ];
      xm[2] = x[idx[m]+2*N];
      phi_m = 0;
      for(int n = 0; n<N; n++)
	{
	  z   = xm[2]-x[n+2*N];
	  rho2= ( (xm[0]-x[n  ])*(xm[0]-x[n  ])+
	          (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  b   = rho2*xi2;
	  qn  = q[n];
	  for(int j2 = -opt.layers; j2<=opt.layers; j2++)
	    {
	      if(j2 == 0)
		continue;
	      
	      k3  = TwoPiOverL*j2;
	      k3z = -k3*z;
	      
	      a   = k3*k3/(4.*xi2);
//	      if(exp(-a-1)/(a*b)<1.e-18)
//		  K0 = 0;
//	      else
		K0  = IncompBesselK0_int(a,b,M,X,W);
	      phi_m += qn*cos(k3z)*K0;
	    }
	}
      phi[m] = phi_m/(opt.box[2]);
    }
}


void SE1P_direct_k0(double* restrict phi,
        const int* restrict idx, int nidx,
        const double* restrict x,
        const double* restrict q, int N,
        const ewald_opts opt)
{
    double phi_m, rho2,xm[2];
    const double xi = opt.xi;
    double egamma = 0.57721566490153286061;

    for(int m=0; m<nidx; m++)
    {
        phi_m=0;
        xm[0] = x[idx[m]    ];
        xm[1] = x[idx[m]+N  ];
        for(int n=0; n<N; n++)
        {
	  rho2 = ( (xm[0]-x[n  ])*(xm[0]-x[n  ]) +
		   (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  if(rho2>__DBL_EPSILON__)
	  	phi_m += -q[n]*(gsl_sf_expint_E1(rho2*xi*xi)+log(rho2*xi*xi)+egamma);
        }
        phi[m] = phi_m/opt.box[2];
    }
}

void SE1P_direct_self(double* restrict phi, 
		      const int* restrict idx, int nidx,
		      const double* restrict q, int N, 
		      const ewald_opts opt)
{
    double c = 2*opt.xi/sqrt(PI);
    for(int m=0; m<nidx; m++)
	phi[m] = -c*q[idx[m]];
} 

void SE1P_direct_force(double* restrict force, 
		       const int* restrict idx, int nidx,
		       const double* restrict x, 
		       const double* restrict q, int N,
		       const ewald_opts opt)
{
  double f[3];
#ifdef _OPENMP
#pragma omp parallel for private(f)
#endif
  for(int m=0; m<nidx; m++)
    {
      f[0] = 0; f[1] = 0; f[2] = 0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p2 = -opt.layers; p2<=opt.layers; p2++)
	    {
	      if(idx[m] == n && p2 == 0)
		continue;

	      double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      double r3 = r*r*r;
	      
	      double c = qn/r3;
	      f[0] += c*rvp[0];
	      f[1] += c*rvp[1];
	      f[2] += c*rvp[2];
	    }
	}
      force[idx[m]    ] = -.5*q[idx[m]]*f[0];
      force[idx[m]+  N] = -.5*q[idx[m]]*f[1];
      force[idx[m]+2*N] = -.5*q[idx[m]]*f[2];
    }
}


void SE1P_direct_real_force(double* restrict force, 
			    const int* restrict idx, int nidx,
			    const double* restrict x, 
			    const double* restrict q, int N,
			    const ewald_opts opt)
{
  double f[3];
  double xi = opt.xi;
  double xi2 = xi*xi;
#ifdef _OPENMP
#pragma omp parallel for private(f)
#endif
  for(int m=0; m<nidx; m++)
    {
      f[0] = 0; f[1] = 0; f[2] = 0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p2 = -opt.layers; p2<=opt.layers; p2++)
	    {
	      if(idx[m] == n && p2 == 0)
		continue;

	      double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      double r2 = r*r;
	      
	      double c = qn*(2*xi/sqrt(PI)*exp(-xi2*r2)+ erfc(xi*r)/r)/r2;
	      f[0] += c*rvp[0];
	      f[1] += c*rvp[1];
	      f[2] += c*rvp[2];
	    }
	}
      force[idx[m]    ] = -.5*q[idx[m]]*f[0];
      force[idx[m]+  N] = -.5*q[idx[m]]*f[1];
      force[idx[m]+2*N] = -.5*q[idx[m]]*f[2];
    }
}

void SE1P_direct_real_rc_force(double* restrict force, 
			       const int* restrict idx, int nidx,
			       const double* restrict x, 
			       const double* restrict q, int N,
			       const ewald_opts opt)
{
  double f[3];
  double xi = opt.xi;
  double xi2 = xi*xi;
#ifdef _OPENMP
#pragma omp parallel for private(f)
#endif
  for(int m=0; m<nidx; m++)
    {
      f[0] = 0; f[1] = 0; f[2] = 0; 
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
      for(int n=0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n], xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double qn = q[n];

	  for(int p2 = -opt.layers; p2<=opt.layers; p2++)
	    {
	      if(idx[m] == n && p2 == 0)
		continue;

	      double rvp[] = {rvec[0], rvec[1],rvec[2]+p2*opt.box[2]};
	      double r = sqrt(rvp[0]*rvp[0]+rvp[1]*rvp[1]+rvp[2]*rvp[2]);
	      double r2 = r*r;

	      if(r > opt.rc) continue;
	      
	      double c = qn*(2*xi/sqrt(PI)*exp(-xi2*r2)+ erfc(xi*r)/r)/r2;
	      f[0] += c*rvp[0];
	      f[1] += c*rvp[1];
	      f[2] += c*rvp[2];
	    }
	}
      force[idx[m]    ] = -.5*q[idx[m]]*f[0];
      force[idx[m]+  N] = -.5*q[idx[m]]*f[1];
      force[idx[m]+2*N] = -.5*q[idx[m]]*f[2];
    }
}

void
SE1P_direct_k0_force(double* restrict force,
		     const int* restrict idx, int nidx,
		     const double* restrict x,
		     const double* restrict q, int N,
		     const ewald_opts opt)
{
  double rho2,xm[2];
  const double xi = opt.xi;
    
  for(int m=0; m<nidx; m++)
    {
      double force_x=0, force_y=0;
      xm[0] = x[idx[m]    ];
      xm[1] = x[idx[m]+N  ];
      for(int n=0; n<N; n++)
        {
	  rho2 = ( (xm[0]-x[n  ])*(xm[0]-x[n  ]) +
		   (xm[1]-x[n+N])*(xm[1]-x[n+N]) );
	  if(rho2==0)
 	     continue;
	  force_x += q[n]*2.*(xm[0]-x[n]  )/rho2*(1-exp(-rho2*xi*xi));
	  force_y += q[n]*2.*(xm[1]-x[n+N])/rho2*(1-exp(-rho2*xi*xi));
        }
      force[idx[m]    ] = -.5*q[idx[m]]*force_x/opt.box[2];
      force[idx[m]+  N] = -.5*q[idx[m]]*force_y/opt.box[2];
      force[idx[m]+2*N] = 0;
    }  
}


void SE1P_direct_fd_force(double* restrict force, 
			  const int* restrict idx, int nidx,
			  const double* restrict x, 
			  const double* restrict q, int N,
			  const ewald_opts opt)
{
  double f[3];
  const double xi   = opt.xi;
  double xi2        = xi*xi;
  double TwoPiOverL = 2.*PI/opt.box[2];
  /* get ready for gsl integration*/
   gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (STACK_SIZE);

  for(int m=0; m<nidx; m++)
    {
      double xm[] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};

      f[0] = 0; f[1] = 0; f[2] = 0;
      for(int n = 0; n<N; n++)
	{
	  double rvec[] = {xm[0]-x[n],xm[1]-x[n+N], xm[2]-x[n+2*N]};
	  double rho2 = rvec[0]*rvec[0] + rvec[1]*rvec[1];
	  double b   = rho2*xi2;
	  double qn  = q[n];
	  for(int j2 = -opt.layers; j2<=opt.layers; j2++)
	    {
	      if(j2 == 0)
		continue;
	      
	      double k3  = TwoPiOverL*j2;
	      double k3z = -k3*rvec[2];
	      
	      double a   = k3*k3/(4.*xi2);
	      //K0  = IncompBesselK0_Simpson(1e-15, &rep, a,b,1);
	      double K0  = call_gsl_bessel_integrator(a,b,w,1);
	      f[0] += 2.*qn*xi2*cos(k3z)*rvec[0]*K0;
	      f[1] += 2.*qn*xi2*cos(k3z)*rvec[1]*K0;
	      //K0  = IncompBesselK0_Simpson(1e-15, &rep, a,b,0);
	      K0  = call_gsl_bessel_integrator(a,b,w,0);
	      f[2] += -qn*k3*sin(k3z)*K0;
	    }
	}
      force[idx[m]    ] = -0.5*q[idx[m]]*f[0]/(opt.box[2]);
      force[idx[m]+  N] = -0.5*q[idx[m]]*f[1]/(opt.box[2]);
      force[idx[m]+2*N] = -0.5*q[idx[m]]*f[2]/(opt.box[2]);
    }

    gsl_integration_workspace_free (w);
}

