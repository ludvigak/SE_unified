#include "stresslet_real_rc.h"

#ifdef BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

#define SWAP(x,y) { tmp=x;x=y;y=tmp; }
static void quicksort(int* restrict list, int* restrict slave, int m, int n);
static void build_cell_list(
			    // Input
			    const double* restrict x, 
			    const int N,
			    const double* restrict box,
			    const double rc,
			    // Output
			    double* rn_p,
			    int ncell[3],
			    int* restrict *ll_p,
			    int* restrict *head_p
			    );


// ==== GENERATE TRIPLETS FOR MATRIX ASSEMBLY
void  get_rs_triplets (const double* restrict x, const double* restrict nvec, const int N,
		       const double* restrict box, const double xi, const double rc, const int nlhs,
		       int* restrict *row_p, int* restrict *col_p, double* restrict val[3][3],
		       int* restrict *buck_size_p, int* restrict *idx_in_array_p, int* numel_p
		       )
{
    // Setup output variables
    int* restrict row;
    int* restrict col;
    int* restrict idx_in_array;
    int* restrict buck_size;

    // Setup variables
    int i,j;
    int ncell[3];
    int* restrict ll;
    int* restrict head;
    double rn;


    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    struct timeval tic, toc;
    gettimeofday(&tic, NULL);
    double time_spent;

    // Build cell list
    build_cell_list(x, N, box, rc, &rn, ncell, &ll, &head);

    if(VERBOSE)
    {
	__PRINTF("[RSRC] SPARSE MATRIX\n");
	__PRINTF("[RSRC] %s, xi=%g\n", OP_TAG, xi);
	__PRINTF("[RSRC] rc=%.3f, rn=%.3f\n", rc, rn);
	__PRINTF("[RSRC] box=(%g,%g,%g), ncell=(%d,%d,%d)\n", 
		  box[0],box[1],box[2],
		  ncell[0],ncell[1],ncell[2]);
    }

    //============================================================
    // CALCULATE INTERACTIONS
    //
    // For all vectors, go through neighbors and save interactions
    // in vectors that are used to create a sparse matrix

    // Allocate a guess based on average density +50%
    int maxel = round(  1.5 * N*N*4*PI*rc*rc*rc/3/(box[0]*box[1]*box[2]) );
    int numel = 0;
    row = __MALLOC(maxel*sizeof(int));
    col = __MALLOC(maxel*sizeof(int));
    for(i=0;i<=2;i++)
	for(j=i;j<=2;j++)
	{
	    val[i][j] = __MALLOC(maxel*sizeof(double));
	}

#ifdef _OPENMP
    int barrier_variable=0;
    int barrier_passed=0;
    int realloc_done=0;
#pragma omp parallel private(i,j) shared(numel,maxel,row,col,val,box,x,nvec,head,ll,px,py,pz,ncell,rn,barrier_variable,barrier_passed,realloc_done) default(none)
#endif
    { // Begin parallel section
    int head_idx;
    int icell[3], home_cell[3];

    int idx_s,idx_t,ip;
    double rsq;
    double pshift[3], xs[3], ns[3], nt[3], xr[3];
    double A1[3][3], A2[3][3];

    const double rcsq = rc*rc;

    // Allocate a bufffer of interactions to be written
    // into triplet list
    const int buf_size = 256;
    int buf_cnt = 0;
    int idx_buf, next_idx_t;
    int* restrict buf_idx_t;
    double* restrict buf_xr;
    double* restrict buf_rsq;
    double* restrict C;
    double* restrict D;

    int tnum = 0;
#ifdef _OPENMP
    tnum = omp_get_thread_num();
    // Seems mxMalloc/mxFree are not thread safe
#pragma omp critical
    {
#endif
    buf_idx_t = __MALLOC(buf_size*sizeof(int));
    buf_xr  = __MALLOC(3*buf_size*sizeof(double));
    buf_rsq = __MALLOC(buf_size*sizeof(double));
    C = __MALLOC(buf_size*sizeof(double));
    D = __MALLOC(buf_size*sizeof(double));

#ifdef _OPENMP
    }
#pragma omp for schedule(dynamic) nowait
#endif
    // Loop over all points
    for(idx_s=0;idx_s<N;idx_s++)
    {
	// Source point
	xs[0] = x[idx_s    ];
	xs[1] = x[idx_s+N  ];
	xs[2] = x[idx_s+2*N];
	// Source point normal vector
	ns[0] = nvec[idx_s    ];
	ns[1] = nvec[idx_s+N  ];
	ns[2] = nvec[idx_s+2*N];
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
	    head_idx = 
		icell[0] +
		icell[1]*ncell[0] + 
		icell[2]*ncell[1]*ncell[0];	    
	    // Go through cell list
	    idx_t = head[head_idx];
	    while(1)
	    {
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
		// Save location of next point in cell chain
		if(idx_t == -1) 
		    next_idx_t = -1;
		else
		    next_idx_t = ll[idx_t];

		// Empty buffer if last point of last neighbour,
		// or buffer full
		if ( (ip==26 && next_idx_t==-1) || buf_cnt==buf_size)
		{
		    // Check if we have enough space to hold buffer contents
		    int idx_write, can_write;
#ifdef _OPENMP
#pragma omp critical
#endif
		    {
			// Check if buffer holds writing space for me
			if(maxel-numel <= 2*buf_cnt)
			    can_write = 0;
			else
			    can_write = 1;
			// Reserve writing space anyway
			idx_write = numel;
			numel += 2*buf_cnt;
		    }
		    if(can_write==0)
		    {
			int alloc_add = buf_size; // How much to add to allocation
#ifdef _OPENMP
			// Barrier
			// Everybody has to wait here before reallocation
			int num_procs = omp_get_num_threads();
			// Allocate more than a fuller buffer for every thread
			alloc_add = num_procs*buf_size; 
			//#pragma omp critical
			//__PRINTF("[%d] Reached barrier (%d,%d) \n",tnum,barrier_variable,num_procs);
#pragma omp critical
			{
			realloc_done = 0; // Everybody agrees reallocation has not been done
			barrier_variable++; // Announce you arrived at barrier
			}
			while(barrier_variable != num_procs) {
#pragma omp flush (barrier_variable) // Spin
			};
#pragma omp atomic
			barrier_passed++; // Anounce you passed barrier
			//#pragma omp critical
			//__PRINTF("[%d] Passed barrier (%d,%d) \n",tnum,barrier_variable,num_procs);
#pragma omp critical			
			{ // Critical section
			if(realloc_done==0) 
			{
			realloc_done=1;
#endif
			// Allocate for full buffer(s) + 20% more
			int new_maxel = ceil(1.2*(maxel+alloc_add));
			__PRINTF("[%d] Reallocating triplet vectors %d -> %d\n",tnum,maxel,new_maxel);
			maxel = new_maxel;
			row = __REALLOC(row, maxel*sizeof(int));
			col = __REALLOC(col, maxel*sizeof(int));
			for(i=0;i<=2;i++)
			    for(j=i;j<=2;j++)
				val[i][j] = __REALLOC(val[i][j], maxel*sizeof(double));
#ifdef _OPENMP
			//__PRINTF("[%d] Done \n",tnum);
			}
			else
			{
			    //__PRINTF("[%d] Someone else reallocated \n",tnum);
			}
			if(barrier_passed==num_procs)
			{
			    //__PRINTF("[%d] Everybody passed! \n",tnum);
			    barrier_passed=0;
			    barrier_variable=0;
			}			
			} // End critical section
#endif
		    } // End can_write==0

		    // Do delayed calculations
		    op_A_CD(C,D,buf_rsq,buf_cnt,xi);

		    //#pragma omp critical
		    //__PRINTF("[%d] Begin write \n",tnum);

		    // Write triplets
		    for(idx_buf=0;idx_buf<buf_cnt;idx_buf++)
		    {
			idx_t = buf_idx_t[idx_buf];
			for(i=0;i<3;i++)
			    xr[i] = buf_xr[3*idx_buf+i];
		
			// Source point normal vector
			nt[0] = nvec[idx_t    ];
			nt[1] = nvec[idx_t+N  ];
			nt[2] = nvec[idx_t+2*N];

			// Calculate interactions t->s and s<-t
			op_A_symm_CD(A1,A2,xr,ns,nt,xi,C[idx_buf],D[idx_buf]);

			// Append results to row,col,val vectors
			row[idx_write] = idx_t;
			col[idx_write] = idx_s;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][idx_write] = A1[i][j];
			    }
			idx_write++;

			row[idx_write] = idx_s;
			col[idx_write] = idx_t;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][idx_write] = A2[i][j];
			    }
			idx_write++;
		    } // endfor buffer

		    //#pragma omp critical
		    //__PRINTF("[%d] End write \n",tnum);

		    buf_cnt = 0;
		} // endif chainend or buffull
		idx_t = next_idx_t;
		if(idx_t == -1)
		    break; // Chain ended
	    } // End of neighbours in this cell
	} // End of cells
    } // End of particles

    //#pragma omp critical
    //__PRINTF("[%d] Exit loop \n",tnum);

#ifdef _OPENMP
#pragma omp atomic
    // If anyone is stuck in barrier, let them know we made it out
    barrier_variable++;
#pragma omp atomic
    barrier_passed++;
#pragma omp critical
#endif
    {
    __FREE(buf_idx_t);
    __FREE(buf_xr);
    __FREE(buf_rsq);
    __FREE(C);
    __FREE(D);
    }
    } // End parallel section

    __FREE(head);
    __FREE(ll);

    gettimeofday(&toc, NULL);
    time_spent = DELTA(tic,toc);

    if(VERBOSE)
    {
	__PRINTF("[RSRC] Triplets generated in %.3f seconds.\n", time_spent);
    }

    //============================================
    // SORT RESULTS WITH COUNTING + QUICK SORT
    // Counting sort on columns, then quicksort on rows
    // in each column
    // (Turns out this is counting sort rather than bucket sort,
    // which I initially thought, hence the buck_* naming.)
    gettimeofday(&tic, NULL);
    buck_size = __MALLOC(N*sizeof(int));
    idx_in_array = __MALLOC(numel*sizeof(int));
    int* restrict buck_count = __MALLOC(N*sizeof(int));
    int* restrict buck_pos = __MALLOC(N*sizeof(int));
    int buck_idx,new_idx;

    // Init lists
    for(i=0;i<N;i++)
    {
	buck_size[i]=0;
	buck_count[i]=0;
    }

    // Count number of elements in each bucket (column)
    for(i=0;i<numel;i++)
    {
	buck_idx = col[i];
	buck_size[buck_idx]++;
    }

    // Cumulative addition to get locations of each bucket after sort,
    // + save largest bucket size for later.
    buck_pos[0] = 0;
    for(i=1;i<N;i++)
    {
	buck_pos[i] = buck_pos[i-1]+buck_size[i-1];
    }

    // Assign each element to a bucket, store permutations in idx_in_array
    int* restrict rowtmp = __MALLOC(numel*sizeof(int));
    for(i=0;i<numel;i++)
    { 
	buck_idx = col[i];
	new_idx = buck_pos[buck_idx] + buck_count[buck_idx];
	idx_in_array[ new_idx ] = i;
	buck_count[buck_idx]++;
    }

    __FREE(buck_count); // Free counter

    // Sort rows using permutations
    // (work-shared)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(i=0;i<numel;i++)
	rowtmp[i] = row[ idx_in_array[i] ];
    __FREE(row);
    row = rowtmp;

    if(nlhs==1)
    {
	__FREE(col); // Free column list if only returning matrix,
    }
    else
    {
	// else sort columns too.
	// Could be done faster with bucket info, but sorted columns are 
	// not needed for real application.
	int* restrict coltmp = __MALLOC(numel*sizeof(int));
	for(i=0;i<numel;i++)
	    coltmp[i] = col[ idx_in_array[i] ];
	__FREE(col);
	col = coltmp;
    }
 
   gettimeofday(&toc,NULL);
   time_spent = DELTA(tic,toc);
   if(VERBOSE)
       __PRINTF("[RSRC] Counting sort of cols finished in %.3f seconds.\n", time_spent);

   gettimeofday(&tic,NULL);

   // Quicksort on buckets
   // Each bucket contains a compressed column.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) shared(buck_pos,buck_size,idx_in_array,row)
#endif
   for(buck_idx=0;buck_idx<N;buck_idx++)
   {
       int begin = buck_pos[buck_idx];
       int size = buck_size[buck_idx];	    
       quicksort(row, idx_in_array, begin, begin+size-1)  ;
   }   
   __FREE(buck_pos); // Free bucket list    

   gettimeofday(&toc,NULL);
   time_spent = DELTA(tic,toc);
   if(VERBOSE)
       __PRINTF("[RSRC] Quicksort of rows finished in %.3f seconds.\n", time_spent);

    // Set return pointers
    *row_p = row;
    *col_p = col;
    *buck_size_p = buck_size;
    *idx_in_array_p = idx_in_array;
    *numel_p = numel;
}

// ==== BUILD CELL LIST
//
// TODO: Add some assertions to make sure rc not too big,
// and that box can be divided into square cells.
static void build_cell_list(
			    // Input
			    const double* restrict x, 
			    const int N,
			    const double* restrict box,
			    const double rc,
			    // Output
			    double* rn_p,
			    int ncell[3],
			    int* restrict *ll_p,
			    int* restrict *head_p
			    )
{
    int i,j;
    int head_idx, ncell_tot;
    int icell[3];
    int* restrict ll;
    int* restrict head;
    double boxmin, rn;

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

    // Prepare cell list
    ll = __MALLOC(N*sizeof(int));
    head = __MALLOC(ncell_tot*sizeof(int));
    for(i=0; i<ncell_tot; i++)
    	head[i] = -1;
    // Do cell partitioning 
    for(i=0; i<N; i++)
    {
	for(j=0; j<3; j++)
	    icell[j] = floor( x[i+N*j]/rn );
	head_idx = 
	    icell[0] +
	    icell[1]*ncell[0] + 
	    icell[2]*ncell[1]*ncell[0];
	ll[i] = head[head_idx];
	head[head_idx] = i;
    }

    *rn_p = rn;
    *ll_p = ll;
    *head_p = head;
}

//============ QUICKSORT ROUTINE
// Applies quicksort on an interval (m,n) of *list, 
// performs the same permutations on *slave.
// Uses a private stack instead of making recursive calls.
static void quicksort(int* restrict list, int* restrict slave, int m, int n) {
    #define  MAX_LEVELS  64
    int beg[MAX_LEVELS], end[MAX_LEVELS]; // Stack
    int key,i,j,k,s,tmp;
    s=0;
    beg[0]=m; 
    end[0]=n;
    while (s>=0) 
    { // While work in stack, pop
	m=beg[s];
	n=end[s];
	if (m<n)
	{
	    k = m+(n-m)/2; // Choose middle for pivot
	    SWAP(list[m],list[k]); // Swap out pivot
	    SWAP(slave[m],slave[k]);
	    // Do quicksort
	    key = list[m];
	    i = m+1;
	    j = n;
	    while(i <= j)
	    {
		while((i <= n) && (list[i] <= key))
		    i++;
		while((j >= m) && (list[j] > key))
		    j--;
		if( i < j)
		{
		    SWAP(list[i],list[j]);
		    SWAP(slave[i],slave[j]);
		}
	    }
	    // Swap in pivot at right place
	    SWAP(list[m],list[j]);
	    SWAP(slave[m],slave[j]);

	    if(s == MAX_LEVELS-1) // Stack full
	    {
		__PRINTF("ERROR. Quicksort reached MAX_LEVELS\n");
		return;
	    }

	    // Recursively sort the lesser list
	    beg[s] = m;
	    end[s] = j-1;
	    beg[s+1]=j+1;
	    end[s+1]=n;
	    s += 1;
	    // Do shortest interval first to limit stack use
	    if (end[s]-beg[s]>end[s-1]-beg[s-1])
	    {
		SWAP(beg[s],beg[s-1]);
		SWAP(end[s],end[s-1]);
	    }
	}
	else 
	{
	    s--; 
	}
    }
}

// ==== Compute result directly
// Do not build sparse matrix
void  compute_rsrc_direct (const double* restrict x, 
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
    int* restrict ll;
    int* restrict head;
    double rn;

    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    struct timeval tic, toc;
    gettimeofday(&tic, NULL);
    double time_spent;

    // Build cell list
    build_cell_list(x, N, box, rc, &rn, ncell, &ll, &head);

    gettimeofday(&toc, NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
    {
	__PRINTF("[RSRC] Cell list built in %.3f seconds.\n", time_spent);
    }

    if(VERBOSE)
    {
	__PRINTF("[RSRC] MATRIX-FREE\n");
	__PRINTF("[RSRC] %s, xi=%g\n", OP_TAG, xi);
	__PRINTF("[RSRC] rc=%.3f, rn=%.3f\n", rc, rn);
	__PRINTF("[RSRC] box=(%g,%g,%g), ncell=(%d,%d,%d)\n", 
		  box[0],box[1],box[2],
		  ncell[0],ncell[1],ncell[2]);
    }

    gettimeofday(&tic, NULL);
#ifdef _OPENMP
#pragma omp parallel shared(phi_out,box,x,nvec,fvec,head,ll,px,py,pz,ncell,rn) default(none)
#endif
    { // Begin parallel section
    // Setup local output
    double* restrict phi = __MALLOC(3*N*sizeof(double));
    for(int i=0;i<3*N;i++)
	phi[i] = 0.0;

    int i,j;
    int head_idx;
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
    // Loop over all points
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
	    head_idx = 
		icell[0] +
		icell[1]*ncell[0] + 
		icell[2]*ncell[1]*ncell[0];	    
	    // Go through cell list
	    idx_t = head[head_idx];
	    while(1)
	    {
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
		// Save location of next point in cell chain
		if(idx_t == -1) 
		    next_idx_t = -1;
		else
		    next_idx_t = ll[idx_t];

		// Empty buffer if last point of last neighbour,
		// or buffer full
		if ( (ip==26 && next_idx_t==-1) || buf_cnt==buf_size)
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
		idx_t = next_idx_t;
		if(idx_t == -1)
		    break; // Chain ended
	    } // End of neighbours in this cell
	} // End of cells
	// Save additions to point s
	phi[idx_s    ] += phi_idx_s[0];
	phi[idx_s+N  ] += phi_idx_s[1];
	phi[idx_s+2*N] += phi_idx_s[2];
    } // End of particles

    // Yes, this reduction is probably crap HPC-wise, 
    // but it works well on my quad core right now.
    for(i=0; i<3*N; i++)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
    	phi_out[i] += phi[i];
    }

#ifdef _OPENMP
    // free/malloc not thread safe under MEX
#pragma omp critical
#endif
    {
    __FREE(phi);
    }

    } // End parallel section

    gettimeofday(&toc, NULL);
    time_spent = DELTA(tic,toc);
    if(VERBOSE)
    {
	__PRINTF("[RSRC] phi computed in %.3f seconds.\n", time_spent);
    }

    *phi_p = phi_out;
}
