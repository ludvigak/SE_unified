#include "../mex/stresslet_real_rc.h"
#include "stresslet_real_rc_new.h"

#ifdef BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

// ******************************** compute_rsrc_direct ******************
// ***********************************************************************
// ==== BUILD CELL LIST (NEW)
//
// TODO: Add some assertions to make sure rc not too big,
// and that box can be divided into square cells.
static void build_cell_list_new(
			    // Input
			    const double* restrict x, 
			    const int N,
			    const double* restrict box,
			    const double rc,
			    // Output
			    double* rn_p,
			    int ncell[3],
			    int* restrict *cell_list_p,
			    int* restrict *cell_idx_p
			    )
{
    int i,j;
    int head_idx, ncell_tot;
    int icell[3];
    double boxmin, rn;
    // Outputs
    int* restrict cell_list;
    int* restrict cell_idx;
    // Intermediates (could do this with fewer vars, but this is clear)
    int* restrict cell_count;
    int* restrict points_in_cell;
    int* restrict point_cell_map;


    // Setup cell partitioning
    boxmin = box[0];

    if(box[1]<boxmin)
	boxmin = box[1];
    if (box[2]<boxmin)
	boxmin = box[2];
    rn = boxmin / floor(boxmin/rc);
    for(i=0;i<3;i++)
	ncell[i] = round( box[i]/rn );
    ncell_tot = ncell[0]*ncell[1]*ncell[2];

    // Prepare arrays
    cell_list  = __MALLOC(N*sizeof(int));
    cell_idx   = __MALLOC((ncell_tot+1)*sizeof(int));
    point_cell_map = __MALLOC(N*sizeof(int));
    points_in_cell = __MALLOC(ncell_tot*sizeof(int));

    for(i=0; i<ncell_tot; i++)
	points_in_cell[i] = 0;

    // Build list in two sweeps 
    for(i=0; i<N; i++)
    {
	for(j=0; j<3; j++)
	    icell[j] = floor( x[i+N*j]/rn );
	int icell_idx = 
	    icell[0] +
	    icell[1]*ncell[0] + 
	    icell[2]*ncell[1]*ncell[0];
	points_in_cell[icell_idx]++;
	point_cell_map[i] = icell_idx;
    }
    // Generate adressing
    cell_idx[0]=0;
    for (int i=0; i<ncell_tot; i++)
	cell_idx[i+1] = cell_idx[i]+points_in_cell[i];
    // Setup new vector
    __FREE(points_in_cell);
    cell_count = __MALLOC(ncell_tot*sizeof(int));
    for(i=0; i<ncell_tot; i++)
	cell_count[i] = 0;
    // Finally build list
    for(i=0; i<N; i++)
    {
	int icell_idx = point_cell_map[i]; 
	int adr = cell_idx[icell_idx] + cell_count[icell_idx];
	cell_list[adr] = i;
	cell_count[icell_idx]++;
    }
    __FREE(cell_count);
    __FREE(point_cell_map);
    *rn_p = rn;
    *cell_list_p = cell_list;
    *cell_idx_p = cell_idx;
}
// ==== Compute result directly
// Do not build sparse matrix
void  compute_rsrc_direct_new (const double* restrict x, 
			   const double* restrict nvec, 
			   const double* restrict fvec, 
			   const int N,
			   const double* restrict box, 
			   const double xi, 
			   const double rc, 
			   double* restrict *phi_p
			   )
{
    // Setup output
    double* restrict phi_out = __MALLOC(3*N*sizeof(double));
    for(int i=0;i<3*N;i++)
	phi_out[i] = 0.0;

    // Setup variables
    int ncell[3];
    int* restrict cell_list;
    int* restrict cell_idx;
    double rn;

    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    struct timeval tic, toc;
    gettimeofday(&tic, NULL);
    double time_spent;

    // Build cell list
    build_cell_list_new(x, N, box, rc, &rn, ncell, &cell_list, &cell_idx);

    gettimeofday(&toc, NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
    {
	__PRINTF("[RSRC] MATRIX-FREE\n");
	__PRINTF("[RSRC] %s, xi=%g\n", OP_TAG, xi);
	__PRINTF("[RSRC] rc=%.3f, rn=%.3f\n", rc, rn);
	__PRINTF("[RSRC] box=(%g,%g,%g), ncell=(%d,%d,%d)\n", 
		  box[0],box[1],box[2],
		  ncell[0],ncell[1],ncell[2]);
	__PRINTF("[RSRC] Cell list built in %.3f seconds.\n", time_spent);
    }
    gettimeofday(&tic, NULL);
#ifdef _OPENMP
#pragma omp parallel \
    shared(phi_out,box,x,nvec,fvec,cell_list,cell_idx,\
	   px,py,pz,ncell,rn) \
    default(none)
#endif
    { // Begin parallel section
    // Setup local output
    double* restrict phi = __MALLOC(3*N*sizeof(double));
    for(int i=0;i<3*N;i++)
	phi[i] = 0.0;

    int i,j;
    int icell_idx;
    int icell[3], home_cell[3];

    int idx_s,idx_t,ip;
    double rsq;
    double pshift[3], xs[3], ns[3], fs[3], nt[3], ft[3], xr[3];
    double A1[3][3], A2[3][3];

    const double rcsq = rc*rc;

 

    // Allocate a bufffer of interactions to be written
    // into triplet list
    const int buf_size = 256;
    int buf_cnt = 0;
    int idx_buf, next_idx_t;
    int buf_idx_t[buf_size];
    double buf_xr[3*buf_size];
    double buf_rsq[buf_size];
    double C[buf_size];
    double D[buf_size];

    int num_procs = 1;
#ifdef _OPENMP
    num_procs = omp_get_num_threads();
    if(VERBOSE)
    {
#pragma omp master	
	__PRINTF("[RSRC] Running on %d threads.\n",num_procs);
    }
#pragma omp for schedule(dynamic) nowait
#endif
    // Loop over all points (work-shared)
    for(idx_s=0;idx_s<N;idx_s++)
    {
	double phi_idx_s[3] = {0.0, 0.0, 0.0};
	// Source point
	xs[0] = x[idx_s    ];
	xs[1] = x[idx_s+N  ];
	xs[2] = x[idx_s+2*N];
	// Source point normal vector
	ns[0] = nvec[idx_s    ];
	ns[1] = nvec[idx_s+N  ];
	ns[2] = nvec[idx_s+2*N];
	// Source point distribution density
	fs[0] = fvec[idx_s    ];
	fs[1] = fvec[idx_s+N  ];
	fs[2] = fvec[idx_s+2*N];
	// Determine home cell
	for(j=0; j<3; j++)
	    home_cell[j] = floor( xs[j]/rn );	
	// Iterate through near cells (including home cell)
	for(ip=0; ip<27; ip++) 
	{
	    // Get neigh cell
	    icell[0] = home_cell[0] + px[ip];
	    icell[1] = home_cell[1] + py[ip];
	    icell[2] = home_cell[2] + pz[ip];
	    // Periodic wrap
	    for(j=0; j<3; j++)
	    {
		// (Could do this with mod)
		pshift[j] = 0;
		if(icell[j] >= ncell[j])
		{
		    icell[j] = 0;
		    pshift[j] = box[j];
		}
		else if(icell[j]<0)
		{
		    icell[j] = ncell[j]-1;
		    pshift[j] = -box[j];
		}
	    }
	    icell_idx = 
		icell[0] +
		icell[1]*ncell[0] + 
		icell[2]*ncell[1]*ncell[0];	    
	    // Go through cell list
	    int cell_a = cell_idx[icell_idx];
	    int cell_b = cell_idx[icell_idx+1];
	    for(int point_idx=cell_a; point_idx<cell_b; point_idx++)
	    {
		idx_t = cell_list[point_idx];
		if(idx_t > idx_s)
		{
		    // r points from s to t
		    for(j=0; j<3; j++)
			xr[j] = x[idx_t+j*N] + pshift[j] - xs[j];
		    // Check if we are within truncation radius
		    rsq = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
		    if(rsq <= rcsq)
		    {		
			// Yes, so put interaction in buffer
			buf_idx_t[buf_cnt] = idx_t;
			buf_rsq[buf_cnt] = rsq;
			for(i=0;i<3;i++)
			    buf_xr[3*buf_cnt+i] = xr[i];
			buf_cnt++;
		    }
		    
		}

		// Empty buffer if last point of last neighbour,
		// or buffer full
		if ( (ip==26 && point_idx==(cell_b-1)) || buf_cnt==buf_size)
		{
		    // Do delayed calculations
		    op_A_CD(C,D,buf_rsq,buf_cnt,xi);

		    // Save interactions
		    for(idx_buf=0;idx_buf<buf_cnt;idx_buf++)
		    {
			idx_t = buf_idx_t[idx_buf];
			for(i=0;i<3;i++)
			    xr[i] = buf_xr[3*idx_buf+i];
		
			// Target point normal vector
			nt[0] = nvec[idx_t    ];
			nt[1] = nvec[idx_t+N  ];
			nt[2] = nvec[idx_t+2*N];

			// Target point distribution density
			ft[0] = fvec[idx_t    ];
			ft[1] = fvec[idx_t+N  ];
			ft[2] = fvec[idx_t+2*N];

			// Calculate interactions t->s and s<-t
			double phi_idx_t[3] = {0.0,0.0,0.0};
			op_A_comp_symm_CD(xr,phi_idx_s,phi_idx_t,ns,nt,fs,ft,xi,C[idx_buf],D[idx_buf]);
			for(i=0; i<3; i++)
			    phi[idx_t+N*i] += phi_idx_t[i];

		    } // endfor buffer

		    buf_cnt = 0;
		} // endif chainend or buffull
	    } // End of neighbours in this cell
	} // End of cells
	// Save additions to point s
	phi[idx_s    ] += phi_idx_s[0];
	phi[idx_s+N  ] += phi_idx_s[1];
	phi[idx_s+2*N] += phi_idx_s[2];
    } // End of particles

#ifdef _OPENMP
    // Yes, this reduction is probably crap HPC-wise, 
    // but it works well on my quad core right now.

    struct timeval tic_red, toc_red;
#pragma omp master
    gettimeofday(&tic_red, NULL);    

    for(i=0; i<3*N; i++)
    {
#pragma omp atomic
    	phi_out[i] += phi[i];
    }

#pragma omp master
    {
	gettimeofday(&toc_red, NULL);    

	double time_spent = DELTA(tic_red,toc_red);
	if(VERBOSE)
	    __PRINTF("[RSRC] Reduction took %.3f seconds.\n", time_spent);
    }

    // free/malloc not thread safe under MEX
#pragma omp critical
    __FREE(phi);
#else
    __FREE(phi_out);
    phi_out = phi;    
#endif

    } // End parallel section

    gettimeofday(&toc, NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
    {
	__PRINTF("[RSRC] phi computed in %.3f seconds.\n", time_spent);
    }

    *phi_p = phi_out;
}


