#include "stresslet_real_rc.h"

#ifdef BEENAKKER
#include "beenakker_op_fd.h"
#else
#error "Must provide -D<method> to compiler"
#endif

// ==== STRUCT AND COMPARISON FUNCTION FOR SORTING
typedef struct
{
    int idx_in_array, row; 
} idx_qsort_buck;

static int compare_qsort_buck(const void* a, const void* b)
{
    int temp = ( * ((idx_qsort_buck*) a) ).row - 
	( * ((idx_qsort_buck*) b) ).row;
    if (temp > 0)
	return 1;
    else if (temp < 0)
	return -1;
    else
	return 0;
}

// ==== GENERATE TRIPLETS FOR ASSEMBLY
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
    int i,j,idx_s,idx_t,ip;
    int head_idx, ncell_tot;
    int ncell[3], icell[3], home_cell[3];
    int* restrict ll;
    int* restrict head;
    double boxmin, rsq;
    double pshift[3], xs[3], ns[3], nt[3], xr[3];
    double A1[3][3], A2[3][3];

    const double rcsq = rc*rc;

    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    clock_t begin, end;
    double time_spent;
    begin = clock();

    //==============================================================
    // BUILD CELL LIST
    //
    // TODO: Add some assertions to make sure rc not too big,
    // and that box can be divided into square cells.

    // Setup cell partitioning
    boxmin = box[0];
    if(box[1]<boxmin)
	boxmin = box[1];
    if (box[2]<boxmin)
	boxmin = box[2];
    const double rn = boxmin / floor(boxmin/rc);
    for(i=0;i<3;i++)
	ncell[i] = round( box[i]/rn );
    ncell_tot = ncell[0]*ncell[1]*ncell[2];

    if(VERBOSE)
    {
	__PRINTF("[RSRC] SPLIT CODE\n");
	__PRINTF("[RSRC] %s, xi=%g\n", OP_TAG, xi);
	__PRINTF("[RSRC] rc=%.3f, rn=%.3f\n", rc, rn);
	__PRINTF("[RSRC] box=(%g,%g,%g), ncell=(%d,%d,%d)\n", 
		  box[0],box[1],box[2],
		  ncell[0],ncell[1],ncell[2]);
    }

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

    // Allocate a bufffer of interactions to be retired later
    const int buf_size = 128;
    int buf_cnt = 0;
    int idx_buf, next_idx_t;
    int* restrict buf_idx_t = __MALLOC(buf_size*sizeof(int));
    double* restrict buf_xr  = __MALLOC(3*buf_size*sizeof(double));
    double* restrict buf_rsq = __MALLOC(buf_size*sizeof(double));
    double* restrict C = __MALLOC(buf_size*sizeof(double));
    double* restrict D = __MALLOC(buf_size*sizeof(double));

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
	    while( idx_t != -1)
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
		next_idx_t = ll[idx_t];

		// Empty buffer if last point of last neighbour,
		// or buffer full
		if ( (ip==26 && next_idx_t==-1) || buf_cnt==buf_size)
		{
		    // Empty buffer by doing delayed calculations
		    op_A_CD(C,D,buf_rsq,buf_cnt,xi);
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

			// Check if we still have allocated space
			if(numel>=maxel-1)
			{
			    __PRINTF("Reallocating...\n");
			    // Else, allocate 20% more
			    maxel = ceil(1.2*maxel);

			    row = __REALLOC(row, maxel*sizeof(int));
			    col = __REALLOC(col, maxel*sizeof(int));
			    for(i=0;i<=2;i++)
				for(j=i;j<=2;j++)
				    val[i][j] = __REALLOC(val[i][j], maxel*sizeof(double));
			}

			row[numel] = idx_t;
			col[numel] = idx_s;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][numel] = A1[i][j];
			    }
			numel++;

			row[numel] = idx_s;
			col[numel] = idx_t;			
			for(i=0; i<=2; i++)
			    for(j=i; j<=2; j++)
			    {
				val[i][j][numel] = A2[i][j];
			    }
			numel++;
		    } // endfor buffer
		    buf_cnt = 0;
		} // endif chainend or buffull
		idx_t = next_idx_t;
	    } // End of neighbours in this cell
	} // End of cells
    } // End of particles

    __FREE(head);
    __FREE(ll);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(VERBOSE)
    {
	__PRINTF("[RSRC] Triplets generated in %.3f seconds.\n", time_spent);
    }

    //============================================
    // SORT RESULTS WITH BUCKET + QUICK SORT
    // Bucket sort on columns, then quicksort on rows
    // in each column
    begin = clock();
    buck_size = __MALLOC(N*sizeof(int));
    idx_in_array = __MALLOC(numel*sizeof(int));
    int* restrict buck_count = __MALLOC(N*sizeof(int));
    int* restrict buck_pos = __MALLOC(N*sizeof(int));
    int buck_idx, buck_max;

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
    buck_max = 0;
    for(i=1;i<N;i++)
    {
	buck_pos[i] = buck_pos[i-1]+buck_size[i-1];
	if(buck_size[i-1]>buck_max)
	    buck_max = buck_size[i-1];
    }

    // Assign each element to a bucket, store permutations in idx_in_array
    for(i=0;i<numel;i++)
    { 
	buck_idx = col[i];
	idx_in_array[ buck_pos[buck_idx] + buck_count[buck_idx] ] = i;
	buck_count[buck_idx]++;
    }

    __FREE(buck_count); // Free counter
    if(nlhs==1)
	__FREE(col); // Free column list if only returning matrix
 
    // Quicksort on buckets (rows indices for a column)
    // Allocate list of pairs to fit largest bucket
    idx_qsort_buck* restrict qlist = __MALLOC(buck_max*sizeof(idx_qsort_buck));
    for(buck_idx=0;buck_idx<N;buck_idx++)
    {
	int begin = buck_pos[buck_idx];
	int size = buck_size[buck_idx];	    
	// Put (idx,row) pairs in list to be qsorted.
	for(i=0; i<size; i++)
	{
	    qlist[i].idx_in_array = idx_in_array[begin+i];
	    qlist[i].row = row[ idx_in_array[begin+i] ];
	}
	qsort(qlist, size, sizeof(idx_qsort_buck), compare_qsort_buck);	    
	for(i=0; i<size; i++)
	    idx_in_array[begin+i] = qlist[i].idx_in_array;
    }       

    __FREE(qlist); // Free temp list
    __FREE(buck_pos); // Free bucket list

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(VERBOSE)
	__PRINTF("[RSRC] Triplets sorted in %.3f seconds.\n", time_spent);

    // Set return pointers
    *row_p = row;
    *col_p = col;
    *buck_size_p = buck_size;
    *idx_in_array_p = idx_in_array;
    *numel_p = numel;
}
