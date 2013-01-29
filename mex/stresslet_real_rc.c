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
			    );
static void barrier(int bar_num, int *barrier_in, int *barrier_out, int *num_procs);
static void transpose(const double* restrict in, double* restrict out, const int N);

// ==== GENERATE TRIPLETS FOR MATRIX ASSEMBLY
void  get_rs_triplets (
		       const double* restrict x_in, 
		       const double* restrict nvec_in, 
		       const int N,
		       const double* restrict box, 
		       const double xi, 
		       const double rc, 
		       const int nlhs,
		       int* restrict *row_p, 
		       int* restrict *col_p, 
		       double* restrict val[3][3],
		       int* restrict *buck_size_p, 
		       int* restrict *idx_in_array_p, 
		       int* numel_p
		       )
{
    // Fix input (legacy format gives bad memory access)
    double* restrict x = __MALLOC(3*N*sizeof(double));
    double* restrict nvec = __MALLOC(3*N*sizeof(double));
    transpose(x_in, x, N);    
    transpose(nvec_in, nvec, N);    

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
    int barrier_in[2]  = {0,0};
    int barrier_out[2] = {0,0};
    int realloc_done=0;
    int num_procs;
#pragma omp parallel private(i,j) \
    shared(numel,maxel,row,col,val,box,x,nvec,head,ll,px,py,pz,ncell,rn,barrier_in,barrier_out,realloc_done,num_procs) \
    default(none)
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
#pragma omp single
    num_procs = omp_get_num_threads();
    if(VERBOSE)
    {
#pragma omp master	
	__PRINTF("[RSRC] Running on %d threads.\n",num_procs);
    }

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
	for(j=0; j<3; j++)
	{
	    // Source point
	    xs[j] = x[idx_s*3+j];
	    // Determine home cell
	    home_cell[j] = xs[j]/rn;	
	    // Source point normal vector
	    ns[j] = nvec[idx_s*3+j];
	}

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
			xr[j] = x[idx_t*3+j] + pshift[j] - xs[j];
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
		    { /* begin critical section */
			// Check if buffer holds writing space for me
			if(maxel-numel <= 2*buf_cnt) {			 
			    can_write = 0;
			    //__PRINTF("[%d] Can't write, reallocation needed! \n",tnum);
			}
			else
			    can_write = 1;
			// Reserve writing in either case
			idx_write = numel;
			numel += 2*buf_cnt;
		    } /* end critical section */

		    /* Begin can_write==0 */
		    if(can_write==0)
		    {
			int alloc_add = buf_size; // How much to add to allocation (single thread)
#ifdef _OPENMP
			// Everybody has to wait here before reallocation
			// Allocate more than a fuller buffer for every thread
			alloc_add = num_procs*buf_size; 
#pragma omp critical
			realloc_done = 0; // Everybody agrees reallocation has not been done
			barrier(0, barrier_in, barrier_out, &num_procs);

#pragma omp critical			
			{ // Critical section			 
			    if(realloc_done==0) 
			    {
				realloc_done=1;
#endif
				// Allocate for full buffer(s) + 20% more
				int new_maxel = ceil(1.2*(maxel+alloc_add));
				if (VERBOSE)
				    __PRINTF("[RSRC][%d] Reallocating triplet vectors %d -> %d\n",tnum,maxel,new_maxel);
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
			}

			barrier(1, barrier_in, barrier_out, &num_procs);
#endif
		    } 
		    /* End can_write==0 */

		    // Do delayed calculations
		    op_A_CD(C,D,buf_rsq,buf_cnt,xi);

		    //#pragma omp critical
		    //__PRINTF("[%d] Begin write \n",tnum);

		    // Write triplets
		    for(idx_buf=0;idx_buf<buf_cnt;idx_buf++)
		    {
			idx_t = buf_idx_t[idx_buf];
			for(i=0;i<3;i++)
			{
			    xr[i] = buf_xr[3*idx_buf+i];
			    // Source point normal vector
			    nt[i] = nvec[idx_t*3+i];
			}

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



#ifdef _OPENMP
#pragma omp critical
    {
	//__PRINTF("[%d] Exit loop , barrier_in={%d,%d}\n",tnum, barrier_in[0], barrier_in[1]);
#pragma omp atomic
    // One less thread going around in loop
    num_procs--;
    }

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
	    icell[j] = floor( x[i*3+j]/rn );
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


//============ Home-brewed barrier
static void barrier(int bar_num, int *barrier_in, int *barrier_out, int *num_procs)
{
#ifdef _OPENMP
    int tnum = omp_get_thread_num();
    // Barrrier arrive
#pragma omp critical
    {
	barrier_in[bar_num]++; // Announce you arrived at barrier
	//__PRINTF("[%d] Reached barrier %d (%d,%d) \n", tnum, bar_num, barrier_in[bar_num], *num_procs);
    }
    // Barrier spin
    while(barrier_in[bar_num] < *num_procs) {
#pragma omp flush
    };
    // Barrier depart
#pragma omp critical			
    {
	barrier_out[bar_num]++; // Anounce you passed barrier
	//__PRINTF("[%d] Passed barrier %d (%d,%d) \n", tnum, bar_num, barrier_out[bar_num], *num_procs);
    }			
    // Barrier reset
#pragma omp critical
    {
	if (barrier_out[bar_num] >= *num_procs)
	{
	    //__PRINTF("[%d] Everybody passed barrier %d. \n",tnum, bar_num);
	    barrier_in[bar_num] = 0;
	    barrier_out[bar_num] = 0;
	}
    }
#endif
}


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
	    icell[j] = x[i*3+j]/rn;
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
// Transpose vector
void transpose(const double* restrict in, double* restrict out, const int N)
{
    for(int i=0; i<N; i++)
    {
	for(int j=0; j<3; j++)
	{
	    out[i*3+j] = in[i+j*N];
	}
    }    
}

// ==== Compute result directly
// Do not build sparse matrix
void  compute_rsrc_direct     (const double* restrict x_in, 
			       const double* restrict nvec_in, 
			       const double* restrict fvec_in, 
			       const int N,
			       const double* restrict box, 
			       const double xi, 
			       const double rc, 
			       double* restrict *phi_p
			       )
{
    struct timeval tic, toc;
    gettimeofday(&tic, NULL);
    double time_spent;

    // Fix input (legacy format gives bad memory access)
    double* restrict x = __MALLOC(3*N*sizeof(double));
    double* restrict nvec = __MALLOC(3*N*sizeof(double));
    double* restrict fvec = __MALLOC(3*N*sizeof(double));
    transpose(x_in, x, N);    
    transpose(fvec_in, fvec, N);    
    transpose(nvec_in, nvec, N);    
    gettimeofday(&toc, NULL);
    double time_tr = DELTA(tic,toc);


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

    // Build cell list
    gettimeofday(&tic, NULL);
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
#pragma omp parallel					\
    shared(phi_out,box,x,nvec,fvec,cell_list,cell_idx,	\
	   px,py,pz,ncell,rn)				\
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
	for(i=0; i<3; i++)
	{
	    // Source point
	    xs[i] = x[idx_s*3+i];
	    // Source point normal vector
	    ns[i] = nvec[idx_s*3+i];
	    // Source point distribution density
	    fs[i] = fvec[idx_s*3+i];
	    // Determine home cell
	    home_cell[i] = xs[i]/rn;	
	}

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
			xr[j] = x[idx_t*3+j] + pshift[j] - xs[j];
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
			{
			    xr[i] = buf_xr[3*idx_buf+i];
			    // Target point normal vector
			    nt[i] = nvec[idx_t*3+i];
			    // Target point distribution density
			    ft[i] = fvec[idx_t*3+i];
			}
			// Calculate interactions t->s and s<-t
			double phi_idx_t[3] = {0.0,0.0,0.0};
			op_A_comp_symm_CD(xr,phi_idx_s,phi_idx_t,ns,nt,fs,ft,xi,C[idx_buf],D[idx_buf]);
			for(i=0; i<3; i++)
			    phi[idx_t*3+i] += phi_idx_t[i];

		    } // endfor buffer

		    buf_cnt = 0;
		} // endif chainend or buffull
	    } // End of neighbours in this cell
	} // End of cells
	// Save additions to point s
	for(int i=0; i<3; i++)
	    phi[idx_s*3+i] += phi_idx_s[i];
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

    gettimeofday(&tic, NULL);
    __FREE(x);
    __FREE(nvec);
    __FREE(fvec);
    double* restrict phi_tr = __MALLOC(3*N*sizeof(double));
    for(int i=0; i<N; i++)
    {
	for(int j=0; j<3; j++)
	{
	    phi_tr[i+j*N] = phi_out[i*3+j];
	}
    }    
    __FREE(phi_out);
    gettimeofday(&toc, NULL);
    time_tr += DELTA(tic,toc);
    if(VERBOSE)
    {
	__PRINTF("[RSRC] Transpose time: %.3f seconds.\n", time_tr);
	__PRINTF("[RSRC] phi computed in %.3f seconds.\n", time_spent);
    }


    *phi_p = phi_tr;
}
