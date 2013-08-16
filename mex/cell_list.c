#include "cell_list.h"
#include "math.h"
#include "mex_compat.h"

// =============================================================================
// ==== BUILD LINKED CELL LIST =================================================
//
// TODO: Add some assertions to make sure rc not too big,
// and that box can be divided into square cells.
void build_linked_cell_list(
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
	// ASSERT(head_idx < ncell_tot, 'cell index out of bounds');
	ll[i] = head[head_idx];
	head[head_idx] = i;
    }

    *rn_p = rn;
    *ll_p = ll;
    *head_p = head;
}

// =============================================================================
// ==== BUILD CELL LIST (NEW) ==================================================
//
// Alternative implementation that returns a straight list
// Reads faster, but requires two sweeps for creation
//
// TODO: Add some assertions to make sure rc not too big,
// and that box can be divided into square cells.
void build_cell_list(
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
    int ncell_tot;
    int icell[3];
    double boxmin, rn;
    // Outputs
    int* restrict cell_list;
    int* restrict cell_idx;
    // Intermediates (could do this with fewer vars, but this is clear)
    int* restrict cell_count;
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
    cell_idx   = __CALLOC(ncell_tot+1, sizeof(int));
    point_cell_map = __MALLOC(N*sizeof(int));

    // Build list in two sweeps,
    // First sweep, count number of points in each cell
    for(i=0; i<N; i++)
    {
	for(j=0; j<3; j++)
	    icell[j] = x[i*3+j]/rn;
	int icell_idx = 
	    icell[0] +
	    icell[1]*ncell[0] + 
	    icell[2]*ncell[1]*ncell[0];
	// ASSERT(icell_idx < ncell_tot, 'cell index out of bounds');
	cell_idx[icell_idx + 1]++;
	point_cell_map[i] = icell_idx;
    }
    // Generate adressing
    cell_idx[0]=0;
    for (int i=0; i<ncell_tot; i++)
	cell_idx[i+1] += cell_idx[i];
    // Setup new vector
    cell_count = __CALLOC(ncell_tot, sizeof(int));
    // Second sweep, build list
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
