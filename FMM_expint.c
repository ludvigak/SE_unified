/*------------------------------------------------------------------------
 * FMM_expint.c
 *------------------------------------------------------------------------
 *Fast Multipole code for particles in the plane. At each particle 
 *z(k) the sum q(m)*ein(alpha*(z(m)-z(k))) is computed over 
 *all k, where q(m) are the charges of the particles, and 
 *ein(v) = expint(v^2)+log(v^2)+GammaEulerConst.
 *  
 *The Matlab syntax is
 *
 *[M1,maxlevel] = mex_M1(z,q,tol,levmax,numthreads)
 * 
 *where
 *  z           - Positions of the particles, complex.
 *  q           - The charges of the particles, real or complex.
 *  tol         - The requested accuracy(a bit pessimistic, tol = 1e-13
 *                is sufficient for machine epsilon accuracy)
 *  levmax      - The max number of box subdivisions before we resort to
 *                slow direct evaluation. This is to be able to handle
 *                (albeit slowly) highly non-uniform distributions of 
 *                particles.
 *  numthreads  - The number of computational threads to use. 
 *  M1          - The output sums
 *  maxlevel    - The final depth of the refinement recursion.
 *
 *Some details:
 *
 *SSE3       - This code uses SSE3-instructions in some inner loops, so
 *             a Pentium 4 processor or later is required. For this to 
 *             be efficient, some arrays need to be aligned on a 32 byte
 *             boundary. To accomplish this, gcc's aligned malloc was
 *             rewritten to use Matlab's mxMalloc (the code is in 
 *             mm_mxmalloc.h). 
 *Threads    - The code is multithreaded. To be sure no race conditions 
 *             occur, mutexes are used. In order not to lock up the entire 
 *             target arrays, a number of mutexes are used. Based on some 
 *             quick testing on only two machines, 64 mutexes are used for 
 *             2 threads and 2048 for 4 threads or above. This can most 
 *             likely be tuned.
 *Memory     - It is hard to exactly calculate exactly how much memory
 *             the code uses but an estimate is:
 *                  ints    : 2*n_particles+5*n_boxes+326
 *                  doubles : 2*n_particles+n_terms*(n_terms+2*ntboxes)+
 *                            6*n_threads*n_terms
 *             Where n_particles is the number of particles, n_boxes 
 *             is the total number of boxes at the finest grid, n_terms is 
 *             the number of terms in the multipole expansion, ntboxes is 
 *             the number of boxes containing > n_terms/2 particles at the 
 *             finest grid and n_threads is the number of threads used.
 *Efficiency - This code runs in O(nlogn) time, for n particles, so 
 *             technically it is not an implementation of the fast
 *             multipole method, which is usually O(n). A simple 
 *             adaptive scheme is used, but there will be problems with
 *             efficiency for highly non-uniform distributions of 
 *             particles.
 *
 *By Rikard Ojala, 2011-03-11.
 */
#include "mex_FMM.h"
#include "expint.h"
#include <stdio.h>
#include <stdlib.h>

/*Error messages unique for this function, the rest of the messages are in
 *mex_FMM.h.*/
#define STD_ERROR "Error in mex_M1: "
#define STD_WARN "Warning in mex_M1: "
#define SYNTAX_ERROR "incorrect syntax. \n\n usage : [M1,maxlevel] = mex_M1(z,q,tol,levmax,numthreads)"

/*Most of the working data is bunched up in larger arrays depending on the
 *type of data and where it is used. These are the offsets of the
 *sub-arrays in the main arrays.*/
/*----------------*/
/*realloc_data(ints):*/
#define BOX_OFFSETS_OFFSET         (0)
#define NPARTICLES_IN_BOX_OFFSET   (nside*nside+1)
#define TAYLOR_BOXES_OFFSET        (2*nside*nside+1)
#define TAYLOR_BOX_NBR_OFFSET      (3*nside*nside+1)
#define ASSIGN_TMP_OFFSET          (4*nside*nside+1)
/*int_data:*/
#define IN_BOX_OFFSET              (0)
#define PARTICLE_OFFSETS_OFFSET    (n_particles)
#define ILIST_X_OFFSET             (2*n_particles)
#define ILIST_Y_OFFSET             (2*n_particles+108)
#define INTERACTION_LIST_OFFSET    (2*n_particles+216)
/*double_data:*/
#define M1_C_OFFSET                (0)
#define CC_OFFSET                  (n_particles)
#define DBL_THDATA_OFFSET          (n_particles+n_terms*n_terms+(32-((n_terms*n_terms)&0xf)))
#define MPOLE_C_OFFSET             (0)
#define TAYLOREXP_C_OFFSET         (n_terms*n_terms)
#define TEMP_C_OFFSET              (5*n_terms*n_terms)
#define DBL_THDATA_SIZE            (5*n_terms*n_terms+n_particles+(32-((n_terms*n_terms)&0xf)))


/*The structure that is passed as the input parameter to MpolesWorker
 *and MpolesWorkerSum.*/
typedef struct {
  double *z_x, *z_y, *q;
  double scale;
  double *thread_data;
  double *double_data, *localexp;
  int *realloc_data, *int_data, *cursquare;
  int n_terms, n_particles, ntboxes, nside, taylor_threshold;
  int maxparticles_in_box;
} MpolesWorkerStruct;

/*The structure that is passed as the input parameter to DirectWorker.*/
typedef struct {
  double *z_x, *z_y, *q;
  double scale;
  double *double_data;
  int *int_data, *realloc_data, *ilist_x, *ilist_y, *cursquare;
  int nside, n_particles;
} DirectWorkerStruct;

/*------------------------------------------------------------------------
 *The main multipole driver function. The function spawns threads that: 
 *
 *1) Computes the multipole expansion of the particles in each box
 *2) Traverses the interaction list of these boxes and, depending on
 *   the number of target particles, either evaluates the interactions 
 *   directly at these targets, evaluates the multipole expansion at
 *   these targets, or converts the multipole expansion into a taylor
 *   expansion centered around the center of the target box, and sums these
 *   together.
 *3) After 1 and 2 are finished for all boxes, a new set of threads 
 *   evaluates the local taylor series in each box with enough particles.
 *------------------------------------------------------------------------
 */
void Mpoles(double *z_x, double *z_y, double* q, double scale,
            double* double_data, int* int_data, int* realloc_data, 
            int n_terms, int nside, int taylor_threshold, int n_particles, 
            int num_threads);

/*------------------------------------------------------------------------
 *When we have reached the level where all boxes contain less than n_terms
 *particles, this function computes direct interaction between nearest
 *neighbors as well as box self-interaction.
 *------------------------------------------------------------------------
 */
void Direct(double *z_x, double *z_y, double* q, double scale,
            double* double_data, int* int_data, int* realloc_data, 
            int nside, int n_particles, int num_threads);

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *z_x, double *z_y, int n_particles, int nside,
            int *maxparticles_in_box, int *int_data, int *realloc_data);

/*------------------------------------------------------------------------
 *This function calls the FMM for evaluting ein sum.
 *------------------------------------------------------------------------
 */
void Do_FMM_Ein(double* M1, int n_particles, double scale,
		double* z_x, double* z_y, double* q,double tol, int numthreads, int lev_max);


#define OUTPUT_OFFSET  (0)
#define ZX_OFFSET      (n_particles)
#define ZY_OFFSET      (2*n_particles)
#define Q_OFFSET       (3*n_particles)
#define ALL_DATA_SZ    (4*n_particles)

/*-----------------------------------------
 *The main entry function
 *-----------------------------------------*/
int main(int argc, char* argv[]) {
  // Input args
  int k;
  double *q,*z_x,*z_y, *all_input;
  double *M1_c,*M2_c;
  double tol, scale;  /*scale is the ewald parameter*box size*/
  int lev_max,n_particles, numthreads;  
  
  if(argc!=5){
    printf("Error! Input arguments are not provided; N scale tol lev_max.\n");
    exit(EXIT_FAILURE);
  }

  /*Get some constants*/
  n_particles = atoi(argv[1]);
  scale       = atof(argv[2]);
  tol         = atof(argv[3]);
  lev_max     = ceil(pow(n_particles,.25)); // lev_max=6; //This shows a better performance
  lev_max     = atoi(argv[4]);
  numthreads  = 8;
  
  /*Get the positions of the particles.*/
  all_input = (double*) _mm_malloc(ALL_DATA_SZ*sizeof(double),32);
  z_x = &all_input[ZX_OFFSET];
  z_y = &all_input[ZY_OFFSET];
  q   = &all_input[Q_OFFSET];
  
  /* Fill the matrix */
  k = 1;
  //    srand(time(NULL));
  for (int i=0; i<n_particles;i++){
    z_x[i] = (double) rand()/RAND_MAX-0.5;
    z_y[i] = (double) rand()/RAND_MAX-0.5;
    q[i] = k;
    k = -k;
  }
  
  /*Create the output vector for direct*/
  M1_c = &all_input[OUTPUT_OFFSET];
  for (int n=0; n<n_particles; n++)
    M1_c[n] = 0;

  /*Do the FMM ein evaluation*/
  Do_FMM_Ein(M1_c, n_particles, scale, z_x, z_y, q, tol, numthreads, lev_max);
  
  /*Create the output vector for direct*/
  M2_c = (double*) _mm_malloc(n_particles*sizeof(double),32);
  for (int n=0; n<n_particles; n++)
    M2_c[n] = 0;
  
  /* DO the exact */
  if(n_particles<1000){
    for (int n=0; n<n_particles; n++)
      {
	for (int m=0; m<n_particles; m++)
	  {
	    if(m==n)
	      continue;
	    double vx = scale*scale*(z_x[m]-z_x[n])*(z_x[m]-z_x[n]);
	    double vy = scale*scale*(z_y[m]-z_y[n])*(z_y[m]-z_y[n]);
	    double vv = expint_log_euler(vx+vy);
	    M2_c[m] += q[n]*vv;
	  }
      }
    double nrm = 0, ns=0;
    for (int m=0; m<n_particles; m++){
      nrm += (M1_c[m]-M2_c[m])*(M1_c[m]-M2_c[m]);
      ns  += (M1_c[m])*(M1_c[m]);
    }
    printf("Err: %g, sum: %g\n", sqrt(nrm/ns), ns);
  }

  M2_c = NULL;
  all_input = NULL;
  _mm_mxFree(all_input);
  _mm_mxFree(M2_c);

}
    
void Do_FMM_Ein(double* M1, int n_particles, double scale,
		double* z_x, double* z_y, double* q, 
		double tol, 
		int numthreads, int lev_max){
  // Input args
  double *double_data;
  double *M1_c, *C,*CC;
  int *int_data,*realloc_data;
  int *ilist_x,*ilist_y;
  int n_terms,taylor_threshold,csize,i,j,k;
  int current_level,nside,maxparticles_in_box;
  
  TIME.STOP = 0.0;
  
  /*Warn user if lev_max is larger than 16.*/
  if(lev_max > 16) {
    lev_max = 16;
  }
  /*We need at least 3 levels.*/
  if(lev_max < 3) {
    lev_max = 3;
  }
  
  /* start timing */
  double start = gettime();
  
  /*The number of terms required in the multipole expansion to get
   *an accuracy of tol.*/
  n_terms = -(int)ceil(log(tol)/LOG2/4);
  
  if( (n_terms%2) != 2)
    n_terms ++;
  
  /*We set the minimum number of taylro terms to 4 
   *for simplicity of computing the taylor and mpole
   *coeffs in compute_mpole_coeff function.*/
  if(n_terms<4)
    n_terms = 4;
  
  /*The breakpoint in terms of number of particles in a box deciding
   *if to use summed up local taylor expansions or direct evaluated
   *multipole expansions.*/
  taylor_threshold = n_terms*4;
  
  /*The main double data vector. It is aligned on a 32 byte boundary
   *to allow for SSE instructions. Contains:
   * M1_c[n_particles]               Output field. Real. Read/Write. Mutexed
   * CC[n_terms*n_terms]             Binomial matrix. Real. Read-only
   * --- One for each thread --- 
   * mpole_c[n_terms*n_terms]        Temp data. Separated from other threads
   * taylorexp_c[4*n_terms*n_terms]  Temp data. Separated from other threads
   * temp_c[n_particles]             Temp data. Separated from other threads
   */
  int ddsz = n_particles+n_terms*n_terms+(32-((n_terms*n_terms)&0xf));
  int ddth = n_particles+5*n_terms*n_terms+(32-((n_terms*n_terms)&0xf));
  double_data = _mm_mxMalloc((ddsz + numthreads*ddth)*sizeof(double),32);
  M1_c = &double_data[M1_C_OFFSET];
  CC = &double_data[CC_OFFSET];
  
  /*Clear the output*/
  for(j=0;j<n_particles;j++)
    M1_c[j] = 0;
  
  /*Compute a matrix containing binomial numbers. It is used when
   *converting from multipole to taylor expansions.*/
  csize = 2*n_terms-1;
  C = calloc(csize*csize,sizeof(double));
  
  for(i = 0;i<csize;i++)
    C[i] = 1;
  
  for(i = 2;i<=csize;i++)
    for(j=2;j<=i;j++)
      C[(j-1)*csize +i-1] = C[(j-2)*csize +i-2]+C[(j-1)*csize +i-2];
  k = 1;
  for(i = 1;i<=n_terms;i++) {
    for (j = 1;j<=n_terms;j++)
      CC[(j-1)*n_terms+i-1]=k*C[(i-1)*csize+j+i-2];
    k*=-1;
  }
  free(C);
  
  /*We start with a level-2 grid, that is 4x4 boxes. On coarser grids
   *the interaction list is empty, so a bit of work would be wasted
   *setting up for these cases.*/
  current_level = 2;
  nside = 4;
  
  /*int_data. Contains data which is read-only in the threads. Contains:
   *
   *in_box[n_particles]           The box point j is assigned to
   *particle_offsets[n_particles] The addresses of the particles in each 
   *                              box.
   *ilist_x[108]                  x-coordinates of the boxes in the 
   *                              interaction list.
   *ilist_y[108]                  y-coordinates of the boxes in the 
   *                              interaction list.
   *interactionlist[108]          interaction list in terms of relative
   *                              box numbers.
   */
  int_data = calloc((2*n_particles+3*108),sizeof(int));
  
  ilist_x = &int_data[ILIST_X_OFFSET];
  ilist_y = &int_data[ILIST_Y_OFFSET];
  /*Initialize the interaction list. See mex_FMM.h.*/
  Intlist(ilist_x,ilist_y);
  
  /*Set up the mutexes. output_mutex[] and localexp_mutex[] are mutexes 
   *protecting the output and localexp_c arrays.*/
  if(numthreads <= 2)
    num_mutexes = 64;
  else
    num_mutexes = 4*512;
  
  output_mutex = malloc(num_mutexes*sizeof(MUTEX_TYPE));
  localexp_mutex = malloc(num_mutexes*sizeof(MUTEX_TYPE));
  for(i=0;i<num_mutexes;i++) {
    INIT_MUTEX(&output_mutex[i]);
    INIT_MUTEX(&localexp_mutex[i]);
  }
  INIT_MUTEX(&mpoleSqTblmutex);
  INIT_MUTEX(&directSqTblmutex);
  
  /*The realloc_data array contains read-only data shared by the 
   *threads. The array is realloc'd each level in the hierarchy since
   *it contains members that depend on the number of sides of the current
   *grid. Calloc-ing this once and for all would be bad since for 
   *uneven distributions of particles many layers may be needed to 
   *resolve tight clusters of particles. A safe size would thus be 
   *prohibitively large. Instead we use a more moderate initial size and
   *start reallocing when this memory is no longer sufficient. 
   *
   *The array contains: 
   *box_offsets[nside*nside+1]          The addresses of the boxes in 
   *                                    particle_offsets
   *nparticles_in_box[nside*nside]      The number of particles in each 
   *                                    box
   *particlenum_in_box[nside*nside+1]   Temporary array in Assign
   *taylor_boxes[nside*nside]           Contains the boxes with more than 
   *                                    taylor_threshold particles. 
   *taylor_box_nbr[nside*nside]         contains the "local numbering" of 
   *                                    the taylor boxes*/
  realloc_data = calloc((5*PREALLOCNSIDE*PREALLOCNSIDE+2),sizeof(int));

  maxparticles_in_box = n_terms+1;
  /*This is the main loop. Work our way down finer grids until the 
   *maximum number of particles in a box is less than a threshold, or  
   *until the maximum number of refinements is reached, in 
   *which case we break and compute the rest of the interactions via  
   *direct evaluation.*/
  while(current_level < lev_max && maxparticles_in_box > n_terms) {
 
    if(nside > PREALLOCNSIDE)
      realloc_data = realloc(realloc_data,(5*nside*nside+2)*sizeof(int));

    Assign(z_x,z_y,n_particles,nside,
	   &maxparticles_in_box,int_data,realloc_data);
    Mpoles(z_x,z_y,q,scale,double_data,int_data,realloc_data,
	   n_terms,nside,taylor_threshold,n_particles,numthreads);

    current_level++;
    nside *= 2;
        
  }
  nside /=2;
  /*Compute the last interactions via direct evaluation.*/
  Direct(z_x,z_y,q,scale,double_data,int_data,realloc_data,nside,n_particles,numthreads);

  /* stop timing */
  double stop = gettime();
  print_log();
    
  /*Clean up mutexes*/
  for(i=0;i<num_mutexes;i++) {
    DESTROY_MUTEX(&output_mutex[i]);
    DESTROY_MUTEX(&localexp_mutex[i]);
  }
  DESTROY_MUTEX(&mpoleSqTblmutex);
  DESTROY_MUTEX(&directSqTblmutex);
    
  /*Free some of the allocated memory*/
  free(output_mutex);
  free(localexp_mutex);
  free(int_data);
  free(realloc_data);
    

  printf("Timing with %d threads: %f\n",numthreads, stop-start);

  /* Copy the result to the output vector */
  for (i=0; i<n_particles; i++)
    M1[i] = M1_c[i];
    
  if(verbose)
    {
      if(DIRECT)
	printf("Direct\n");
      if(TAYLOR)
	printf("Taylor\n");
      if(MPOLE)
	printf("Mpole\n");
    }
    
  /*Free the last of the allocated memory. (Aligned free)*/
  _mm_mxFree(double_data);
  printf("DONE\n");
}

/*------------------------------------------------------------------------
 *The main multipole driver function.
 *The function spawns threads that:
 *1) Computes the multipole expansion of the particles in each box
 *2) Traverses the interaction list of these boxes and, depending on
 *   the number of target particles, either evaluates the interactions 
 *   directly at these targets, evaluates the multipole expansion at
 *   these targets, or converts the multipole expansion into a taylor
 *   expansion centered around the center of the target box, and sums these
 *   together.
 *3) After 1 and 2 are finished for all boxes, a new set of threads 
 *   evaluates the local taylor series in each box with enough particles.
 *------------------------------------------------------------------------
 */
void Mpoles(double *z_x, double *z_y, double* q, double scale,
            double* double_data, int* int_data, int* realloc_data, 
            int n_terms, int nside, int taylor_threshold, int n_particles, 
            int num_threads) {

    MpolesWorkerStruct* arguments;
    /*See below for explanations of these pointers.*/
    int *interaction_list = &int_data[INTERACTION_LIST_OFFSET];
    int *taylor_box_nbr = &realloc_data[TAYLOR_BOX_NBR_OFFSET];
    int *taylor_boxes = &realloc_data[TAYLOR_BOXES_OFFSET];
    double *localexp_c;
    double *double_offset;
    /*The interaction list.*/
    int *ilist_x = &int_data[ILIST_X_OFFSET];
    int *ilist_y = &int_data[ILIST_Y_OFFSET];
    /*The number of particles in each box.*/
    int *nparticles_in_box = &realloc_data[NPARTICLES_IN_BOX_OFFSET];
    /*ntboxes is the number of Taylor boxes and cursquare is the next
     *square to be treated by the threads.*/
    int ntboxes,cursquare=0,i;

    /*Get boxes that contain enough particles to qualify for taylor series
     *conversion. Store the numbering of these in taylor_boxes. ntboxes 
     *is the number of such boxes.*/
    ntboxes=0;
    for(i=0;i<nside*nside;i++)
        if(nparticles_in_box[i] > taylor_threshold)
            taylor_boxes[ntboxes++] = i;

    /*taylor_box_nbr contains the "local numbering" in terms of 
     *taylor_boxes.*/
    for(i=0;i<ntboxes;i++)
        taylor_box_nbr[taylor_boxes[i]] = i;           

    /*localexp_c. Contains the local taylor expansions which are common
     *to all threads.*/
    if(ntboxes > 0)
        localexp_c = _mm_mxMalloc(n_terms*n_terms*ntboxes*sizeof(double),32);
    else
        localexp_c = NULL;

    /*Set up the interaction list.*/
    for(i=0;i<108;i++)  {
        interaction_list[i] = ilist_y[i]*nside+ilist_x[i];
    }
    /*Clear the local Taylor expansions.*/
    for(i=0;i<n_terms*n_terms*ntboxes;i++)
        localexp_c[i] = 0;
    
    /*Alloc some arrays. mpoleWorkerThd holds the thread instances and 
     *arguments holds the call arguments to the threads.*/
    mpoleWorkerThd = malloc(num_threads*sizeof(THREAD_TYPE));
    arguments = malloc(num_threads*sizeof(MpolesWorkerStruct));
    
    /*Keep track of the offsets of the thread's private data.*/
    double_offset = &double_data[DBL_THDATA_OFFSET];
    
    /*Fill the arguments structs and spawn the threads that are 
     *responsible for step 1 and 2 in the description above.*/
    for(i=0;i<num_threads;i++) {       
        arguments[i].z_x              = z_x;
        arguments[i].z_y              = z_y;
        arguments[i].q                = q;
        arguments[i].scale            = scale;
        arguments[i].thread_data      = double_offset;
        arguments[i].double_data      = double_data;
        arguments[i].realloc_data     = realloc_data;
        arguments[i].localexp         = localexp_c;
        arguments[i].int_data         = int_data;
        arguments[i].n_terms          = n_terms;
        arguments[i].n_particles      = n_particles;
        arguments[i].ntboxes          = ntboxes;
        arguments[i].nside            = nside;
        arguments[i].cursquare        = &cursquare;
        arguments[i].taylor_threshold = taylor_threshold;
        
        /*Spawn thread*/
	THREAD_CREATE(mpoleWorkerThd[i],MpolesWorker,(void*) &arguments[i]);
        
        /*Step forward in the private array*/
        double_offset = &double_offset[DBL_THDATA_SIZE];
    }

    /*Wait for all threads to complete */
	for(i = 0;i<num_threads;i++)
        THREAD_JOIN(mpoleWorkerThd[i]);

    /*Start again from the first box*/
    cursquare = 0;
    /*Spawn the threads that do step 3 above.*/  
    for(i=0;i<num_threads;i++)
      THREAD_CREATE(mpoleWorkerThd[i],MpolesWorkerSum,(void*) &arguments[i]);
  
    /*Wait for all threads to complete */
    for(i = 0;i<num_threads;i++)
      THREAD_JOIN(mpoleWorkerThd[i]);
    
    /*Clean up and free the memory used.*/
    _mm_mxFree(localexp_c);
    free(mpoleWorkerThd);
    free(arguments);
  
}
/*------------------------------------------------------------------------
 *Threaded worker function that computes multipole expansions for the
 *particles in each box and then traverses the interaction list doing one
 *of three things: either the interactions are computed directly, or the 
 *multipole expansion is evaluated directly in the target box, or the 
 *multipole expansions are converted into local taylor expansions around 
 *the centers of the target boxes and summed up.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorker(void *argument) {
    /*Get the variables from the argument structure.*/
    MpolesWorkerStruct* arg = (MpolesWorkerStruct*) argument;
    /*The number of terms in the expansions.*/
    int n_terms    = arg->n_terms;
    int twon_terms = 2*n_terms;
    int n2_terms   = n_terms*n_terms;
    /*The total number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    /*The threshold number of particles when we switch to Taylor series.*/
    int taylor_threshold = arg->taylor_threshold;
    
    /*The "global numbering" of the Taylor boxes.*/
    int *taylor_box_nbr = &(arg->realloc_data[TAYLOR_BOX_NBR_OFFSET]);
    /*The relative positions of the interaction list.*/
    int *interaction_list = &(arg->int_data[INTERACTION_LIST_OFFSET]);
    /*The current box to be treated. Common to all threads.*/
    int *cursquare = arg->cursquare;
    
    /*Pointers to particle position and charges.*/
    double *z_x = arg->z_x;
    double *z_y = arg->z_y;
    double *q = arg->q;
    double scale = arg->scale;
    double scale2= scale*scale;
    /*Pointer to the output sums.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    /*A binomial matrix. Used when going from multipole to Taylor.*/
    double *CC = &(arg->double_data[CC_OFFSET]);
    /*Matrix containing the Taylor series for each Taylor box. Common to
     *all threads.*/
    double *localexp_c = arg->localexp;
    /*--- The following arrays are unique to each thread. ---*/
    /*The multipole expansion for the boxes.*/
    double *mpole_c = &(arg->thread_data[MPOLE_C_OFFSET]);
    /*Temporary taylor array.*/
    double *taylorexp_c = &(arg->thread_data[TAYLOREXP_C_OFFSET]);
    /*Buffer array for multipole evaluation.*/
    double *temp_c = &(arg->thread_data[TEMP_C_OFFSET]);
    
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    /*The interaction list.*/
    int *ilist_x = &(arg->int_data[ILIST_X_OFFSET]);
    int *ilist_y = &(arg->int_data[ILIST_Y_OFFSET]);

    double box_center_re,box_center_im;
    /* mpole_v stores the powers
     * the summation over \Sum qi(xi-xA)^k1(yi-yA)^k2.*/
    double* mpole_v = _mm_mxMalloc(n2_terms*sizeof(double),32);
    
    /* Temporary array to store powers of zx and zy .*/
    double* zx_pow = _mm_mxMalloc(n_terms*sizeof(double),32);
    double* zy_pow = _mm_mxMalloc(n_terms*sizeof(double),32);

    /*Loop through the boxes.*/
    for(;;) {
        int j,current_box,i;

        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&mpoleSqTblmutex);
        i = *cursquare;
        
        while(i < nside*nside && nparticles_in_box[i] == 0)
            i++;
        
        /*If we are done we exit immediately.*/
        if(i==nside*nside) {
            UNLOCK_MUTEX(&mpoleSqTblmutex);
            break;
        }

        *cursquare = i+1;
        UNLOCK_MUTEX(&mpoleSqTblmutex);
        
        current_box = i;
        
        /*The center of the current box.*/
        box_center_re = (double)(current_box%nside-(double)(nside-1)/2)/nside;
        box_center_im = (double)(current_box/nside-(double)(nside-1)/2)/nside;

        /*We need to reset the multipole expansion for this box.*/
        for(j=0;j<n2_terms;j++)
            mpole_v[j] = 0;
        
        /*Compute the multipole expansion around the center of this box
         *This loop shows up as a hotspot in Shark, so SIMD instructions
         *are used to speed it up.*/
         for(j=0;j<nparticles_in_box[current_box];j++) {
	   unsigned int k1, k2;
             
            double zxc = box_center_re-z_x[particle_offsets[box_offsets[current_box]+j]];
            double zyc = box_center_im-z_y[particle_offsets[box_offsets[current_box]+j]];
            
            double qj  = q[particle_offsets[box_offsets[current_box]+j]];
            double sx,sy;
	    
            /* compute powers of zx and zy */
	    zxc *= scale;
	    zyc *= scale;
            zx_pow[0]=1.0; zy_pow[0]=1.0;
            for (k1=1; k1<n_terms; k1++){
                zx_pow[k1] = zx_pow[k1-1]*zxc;
                zy_pow[k1] = zy_pow[k1-1]*zyc;
            }

	    /*compute the rest of mpole_a coeffs and add up to form mpole_c.
	     *This loop accounts for 10% of the runtime.*/    
	    for (k1=0; k1<n_terms; k1++){
	      sx = qj*zx_pow[k1];
	      for (k2=0; k2<n_terms-k1; k2++){
		sy = zy_pow[k2];
		mpole_v[k1*n_terms+k2] += sx*sy;
	      }
	    } 
	 }
        
        /*Loop though the interaction list.*/
        for(j=0;j<27;j++) {
            /*Make sure that the target box is inside the computational
             *grid. The ci_flag variable reflects the fact that the inter-
             *action list of adjacent boxes may be similar if they are both
             *in a 2x2 superbox properly aligned. See Beatson and 
             *Greengard page 19. current_box%nside and current_box/nside 
             *are the box coordinates.*/
            int ci_flag = 2*((current_box%nside)&1)+((current_box/nside)&1);
         
            if((current_box%nside + ilist_x[j+ci_flag*27])<nside &&
               (current_box/nside + ilist_y[j+ci_flag*27])<nside &&
               (current_box%nside + ilist_x[j+ci_flag*27])>=0 &&
               (current_box/nside + ilist_y[j+ci_flag*27])>=0) {
                
                /*The number of the target box.*/
                int target_box = current_box+interaction_list[j+ci_flag*27];

                /*Check the number of sources of the target box.*/
                if(nparticles_in_box[target_box] > taylor_threshold) {
 		    TAYLOR = 1;
  		    unsigned int k1, k2, l1, l2;
                    /*If larger than taylor_threshold we convert the
                     *multipole expansion into a local taylor expansion 
                     *centered around the target box center.*/

                    double* tptr2 = &localexp_c[taylor_box_nbr[target_box]*n2_terms];
		    /*The center of the target box. */
		    double tbox_center_re = (double)(target_box%nside-(double)(nside-1)/2)/nside;
		    double tbox_center_im = (double)(target_box/nside-(double)(nside-1)/2)/nside;
                    
                    /*Compute the local taylor expansion. First we compute 
		     *mpole coefficients on a larger set K+L 4n_terms^2 terms
		     *around xA-xB where xA is the source center and xB is the
		     *target center and store it in taylorexp_c. Then with 
		     *multiplication with CC and mpole_v and is summed up with
		     *the previous localexp_c values.*/
		    double t_x = tbox_center_re - box_center_re;
		    double t_y = tbox_center_im - box_center_im;

		    compute_mpole_c(taylorexp_c, twon_terms, scale*t_x, scale*t_y);

		    /*Sum the local taylor expansions; these are common to
                     *all threads, hence the mutex. Actually, there are
                     *several mutexes(tunable), one for each group of
                     *boxes, so that we lock up as little of the localexp-
                     *matrix as possible.*/
                    LOCK_MUTEX(&localexp_mutex[target_box&(num_mutexes-1)]);
		    for (l1=0; l1<n_terms; l1++)
		      for (l2=0; l2<n_terms-l1; l2++){
			double tmp=0;
			  for (k1=0; k1<n_terms; k1++)
			    for (k2=0; k2<n_terms-k1; k2++){
			      double cc1 = CC[k1*n_terms+l1];
			      double cc2 = CC[k2*n_terms+l2];
			      tmp += cc1*cc2*taylorexp_c[(k1+l1)*twon_terms+(k2+l2)]*mpole_v[k1*n_terms+k2];
			    }
			  tptr2[l1*n_terms+l2] += tmp;
		      }
                    UNLOCK_MUTEX(&localexp_mutex[target_box&(num_mutexes-1)]);

                    
                }else if(nparticles_in_box[target_box] > 0) {

		  if(nparticles_in_box[current_box] < 32) {
  		        DIRECT = 1;
		        /*There are relatively few particles in the target
                         *box, and if there are few particles in the 
                         *current box, we may as well evaluate interactions 
                         *directly.*/
                        unsigned int k;
                        int* tptr  = &particle_offsets[box_offsets[current_box]]; 
                        int* tptr2 = &particle_offsets[box_offsets[target_box]]; 
                        LOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                        /*This loop computes the direct interaction between
			 *all the particles in current_box and the 
			 *particles in target_box. This loop accounts for more than
			 *50% of the runtime, so we do SIMD + unrolling.*/

                        for(k=0;k<(unsigned int)nparticles_in_box[target_box];k++) {
			  unsigned int l=0,l0=0;
			  double res=0;
			  double tz_x, tz_y;
			  double cz_x = z_x[tptr2[k]];
			  double cz_y = z_y[tptr2[k]];
#ifdef __AVX__
			  double resv[4] MEM_ALIGN;
			  __m256d ZX,ZY,Ql,V;
			  __m256d SCALE2 = _mm256_set1_pd(scale2);
			  __m256d RESV = _mm256_setzero_pd();
			  __m256d CZX = _mm256_set1_pd(cz_x);
			  __m256d CZY = _mm256_set1_pd(cz_y);

			  for(l=0;l<(unsigned int)nparticles_in_box[current_box]/8*8;l+=8) {
			    ZX = _mm256_set_pd(z_x[tptr[l+3]],z_x[tptr[l+2]],z_x[tptr[l+1]],z_x[tptr[l]]);
			    ZY = _mm256_set_pd(z_y[tptr[l+3]],z_y[tptr[l+2]],z_y[tptr[l+1]],z_y[tptr[l]]);
			    Ql = _mm256_set_pd(q[tptr[l+3]],q[tptr[l+2]],q[tptr[l+1]],q[tptr[l]]);
			    
			    ZX = _mm256_sub_pd(CZX,ZX);
			    ZX = _mm256_mul_pd(ZX,ZX);
			    ZY = _mm256_sub_pd(CZY,ZY);
			    ZY = _mm256_mul_pd(ZY,ZY);
			    V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
			    V = _mm256_mul_pd(Ql,V);
			    RESV = _mm256_add_pd(RESV, V);
			    
			    l0 = 4;
			    ZX = _mm256_set_pd(z_x[tptr[l+3+l0]],z_x[tptr[l+2+l0]],z_x[tptr[l+1+l0]],z_x[tptr[l+l0]]);
			    ZY = _mm256_set_pd(z_y[tptr[l+3+l0]],z_y[tptr[l+2+l0]],z_y[tptr[l+1+l0]],z_y[tptr[l+l0]]);
			    Ql = _mm256_set_pd(q[tptr[l+3+l0]],q[tptr[l+2+l0]],q[tptr[l+1+l0]],q[tptr[l+l0]]);
			    
			    ZX = _mm256_sub_pd(CZX,ZX);
			    ZX = _mm256_mul_pd(ZX,ZX);
			    ZY = _mm256_sub_pd(CZY,ZY);
			    ZY = _mm256_mul_pd(ZY,ZY);
			    V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
			    V = _mm256_mul_pd(Ql,V);
			    RESV = _mm256_add_pd(RESV, V);
			  }
			  l0=l;
			  for(l=l0;l<(unsigned int)nparticles_in_box[current_box]/4*4;l+=4) {
			    ZX = _mm256_set_pd(z_x[tptr[l+3]],z_x[tptr[l+2]],z_x[tptr[l+1]],z_x[tptr[l]]);
			    ZY = _mm256_set_pd(z_y[tptr[l+3]],z_y[tptr[l+2]],z_y[tptr[l+1]],z_y[tptr[l]]);
			    Ql = _mm256_set_pd(q[tptr[l+3]],q[tptr[l+2]],q[tptr[l+1]],q[tptr[l]]);
			    
			    ZX = _mm256_sub_pd(CZX,ZX);
			    ZX = _mm256_mul_pd(ZX,ZX);
			    ZY = _mm256_sub_pd(CZY,ZY);
			    ZY = _mm256_mul_pd(ZY,ZY);
			    V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
			    V = _mm256_mul_pd(Ql,V);
			    RESV = _mm256_add_pd(RESV, V);
			  }
			  l0 = l;
			  for(l=l0;l<(unsigned int)nparticles_in_box[current_box];l++) {
			    tz_x = (cz_x-z_x[tptr[l]])*(cz_x-z_x[tptr[l]]);
			    tz_y = (cz_y-z_y[tptr[l]])*(cz_y-z_y[tptr[l]]);
			    res += q[tptr[l]]*expint_log_euler( scale2*(tz_x+tz_y));
			  }
			  _mm256_store_pd(resv,RESV);
			  res += sum4(resv);
#elif defined __SSE4_2__
			  double resv[2] MEM_ALIGN;
			  __m128d ZX,ZY,Ql, V;
			  __m128d SCALE2 = _mm_set1_pd(scale2);
			  __m128d RESV = _mm_setzero_pd();
			  __m128d CZX = _mm_set1_pd(cz_x);
			  __m128d CZY = _mm_set1_pd(cz_y);
			  for(l=0;l<(unsigned int)nparticles_in_box[current_box]/4*4;l+=4) {
			    ZX = _mm_set_pd(z_x[tptr[l+1]],z_x[tptr[l]]);
			    ZY = _mm_set_pd(z_y[tptr[l+1]],z_y[tptr[l]]);
			    Ql = _mm_set_pd(q[tptr[l+1]],q[tptr[l]]);
			    
			    ZX = _mm_sub_pd(CZX,ZX);
			    ZX = _mm_mul_pd(ZX,ZX);
			    ZY = _mm_sub_pd(CZY,ZY);
			    ZY = _mm_mul_pd(ZY,ZY);
			    V = _mm_ein_pd(_mm256_mul_pd(SCALE2,_mm_add_pd(ZX,ZY)));
			    V = _mm_mul_pd(Ql,V);
			    RESV = _mm_add_pd(RESV, V);
			    
			    l0 = 2;
			    ZX = _mm_set_pd(z_x[tptr[l+1+l0]],z_x[tptr[l+l0]]);
			    ZY = _mm_set_pd(z_y[tptr[l+1+l0]],z_y[tptr[l+l0]]);
			    Ql = _mm_set_pd(q[tptr[l+1+l0]],q[tptr[l+l0]]);
			    
			    ZX = _mm_sub_pd(CZX,ZX);
			    ZX = _mm_mul_pd(ZX,ZX);
			    ZY = _mm_sub_pd(CZY,ZY);
			    ZY = _mm_mul_pd(ZY,ZY);
			    V = _mm_ein_pd(_mm256_mul_pd(SCALE2,_mm_add_pd(ZX,ZY)));
			    V = _mm_mul_pd(Ql,V);
			    RESV = _mm_add_pd(RESV, V);
			  }
			  l0 = l;
			  for(l=l0;l<(unsigned int)nparticles_in_box[current_box];l++) {
			    tz_x = (cz_x-z_x[tptr[l]])*(cz_x-z_x[tptr[l]]);
			    tz_y = (cz_y-z_y[tptr[l]])*(cz_y-z_y[tptr[l]]);
			    res += q[tptr[l]]*expint_log_euler(scale2*(tz_x+tz_y));
			  }
			  _mm_store_pd(resv,RESV);
			  res += resv[0]+resv[1];
#else
			  for(l=0;l<(unsigned int)nparticles_in_box[current_box];l++) {
			    tz_x = (cz_x-z_x[tptr[l]])*(cz_x-z_x[tptr[l]]);
			    tz_y = (cz_y-z_y[tptr[l]])*(cz_y-z_y[tptr[l]]);
			    res += q[tptr[l]]*expint_log_euler(scale2*(tz_x+tz_y));
			  }
#endif
			  M1_c[tptr2[k]] += res;
			}
			UNLOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
		  }else{
		        MPOLE = 1;
                        /*There are too many particles in the current box,
                         *and direct evaluation would be expensive. Instead
                         *evaluate the multipole expansion at each particle
                         *in the target box.*/
			unsigned int k,k0;
                        int* tptr = &particle_offsets[box_offsets[target_box]];
			
                        /*For each particle in the target box compute the multipole
                         *expansion store in mpole_c and using mpole_v evaluate
                         *expansions. For simiplicity mpole_v and mpole_c
			 *have n_terms^2 elements but we only fill those elements
			 *where row+col=n_terms. This should give us two times 
			 *speed-up with the expnse of more memory. (tunable). */
                        for(k=0;k<(unsigned int)nparticles_in_box[target_box];k++) {
			  unsigned int k1,k2;
			  double t_x = (z_x[tptr[k]]-box_center_re);
			  double t_y = (z_y[tptr[k]]-box_center_im);
			  
			  /* compute mpole_c coefficients.*/
			  compute_mpole_c(mpole_c, n_terms, scale*t_x, scale*t_y);
			  temp_c[k] = 0;
			  double tmp = 0;
			  for (k1=0; k1<n_terms; k1++){
			      for (k2=0; k2<n_terms-k1; k2++){
				tmp += mpole_c[k1*n_terms+k2]*mpole_v[k1*n_terms+k2];
			      }
			  }
			  temp_c[k] = tmp;
			}

                        /*Add the contributions, this is safe since we
                         *know that there is less than taylor_threshold
                         *particles in the target box. taylor_threshold is 
                         *n_terms*2 by default.*/
			LOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                        for(k=0;k<(unsigned int)nparticles_in_box[target_box]/4*4;k+=4) {
			  k0=0; M1_c[tptr[k+k0]] += temp_c[k+k0];
			  k0++; M1_c[tptr[k+k0]] += temp_c[k+k0];
			  k0++; M1_c[tptr[k+k0]] += temp_c[k+k0];
			  k0++; M1_c[tptr[k+k0]] += temp_c[k+k0];
                        }
			k0 = k;
                        for(k=k0;k<(unsigned int)nparticles_in_box[target_box];k++) {
			  M1_c[tptr[k]]   += temp_c[k];
                        }
			UNLOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                    }/*Direct or multipole in target box.*/
                } /*Taylor or Direct/Multipole in target box.*/
            }/*Box is inside the computational domain.*/
        }/*Interaction list loop.*/
    }/*Box loop.*/
    _mm_mxFree(mpole_v);
    _mm_mxFree(zx_pow);_mm_mxFree(zy_pow);
    THREAD_EXIT();
}
/*------------------------------------------------------------------------
 *Threaded worker function that evaluates and sums the local taylor 
 *expansions that were computed in MpolesWorker.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorkerSum(void* argument) {
    /*Get the variables from the argument structure.*/
    MpolesWorkerStruct* arg = (MpolesWorkerStruct*) argument;
    /*Number of terms in the expansions.*/
    int n_terms = arg->n_terms;
    /*The total number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    /*Number of Taylor boxes.*/
    int ntboxes = arg->ntboxes;
    
    /*Pointers to particle positions.*/
    double *z_x = arg->z_x,*z_y = arg->z_y;
    double scale = arg->scale;
    /*Pointer to matrix of Taylor series for boxes with enough particles.*/
    double *localexp_c = arg->localexp;
    /*Pointer to the output data.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    
    /*The next box to be treated. Is common to all threads. Has mutex guard.*/
    int *cursquare = arg->cursquare;
    /*The "global numbering" of the boxes with Taylor series.*/
    int *taylor_boxes = &(arg->realloc_data[TAYLOR_BOXES_OFFSET]);
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    
    double box_center_re,box_center_im;
    double *zx_pow = _mm_mxMalloc(n_terms*sizeof(double),32);
    double *zy_pow = _mm_mxMalloc(n_terms*sizeof(double),32);

    /*Loop through the boxes with more than taylor_threshold points in them.*/
    for(;;) {
      int j,current_box,i;
      int* tptr;
      double *tptr2;

        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&mpoleSqTblmutex);
        i = *cursquare;
        /*If we are done we exit immediately.*/
        if(i==ntboxes) {
            UNLOCK_MUTEX(&mpoleSqTblmutex);
            break;
        }

        *cursquare = i+1;
        UNLOCK_MUTEX(&mpoleSqTblmutex);

        /*The number of the current box.*/
        current_box = taylor_boxes[i];
    
        /*The center of the current box.*/
        box_center_re = (double)(current_box%nside-(double)(nside-1)/2)/nside;
        box_center_im = (double)(current_box/nside-(double)(nside-1)/2)/nside;
        
        /*A temporary pointer to the pre-offsets of this box's particles.*/
        tptr = &particle_offsets[box_offsets[current_box]];
	tptr2 = &localexp_c[i*n_terms*n_terms];
	
	
        /*Evaluate the taylor series using Horner's rule.
         *(Minor)Shark hotspot, so we go SIMD.*/
        for(j=0;j<nparticles_in_box[current_box];j++) {
	  unsigned int k1, k2;

	  double zxc = box_center_re - z_x[tptr[j]];
	  double zyc = box_center_im - z_y[tptr[j]];

	  /* compute powers of zx and zy */
	  zxc *= scale;
	  zyc *= scale;
	  zx_pow[0]=1.0; zy_pow[0]=1.0;
	  for (k1=1; k1<n_terms; k1++){
	    zx_pow[k1] = zx_pow[k1-1]*zxc;
	    zy_pow[k1] = zy_pow[k1-1]*zyc;
	  }

	  /*Add the contribution of the taylor expansion to the end 
	   *result. Accesses to M1_c here are completely parallel, so no 
	   *mutexes are necessary.*/
	  for (k1=0; k1<n_terms; k1++)
	    for (k2=0; k2<n_terms-k1; k2++){
	      M1_c[tptr[j]] += zx_pow[k1]*zy_pow[k2]*tptr2[k1*n_terms+k2];
	    }
        }
    }
    _mm_mxFree(zx_pow);
    _mm_mxFree(zy_pow);
    THREAD_EXIT();
    
}

/*------------------------------------------------------------------------
 *When we have reached the level where all boxes contain less than n_terms
 *particles, this function computes direct interaction between nearest
 *neighbors as well as box self-interaction.
 *------------------------------------------------------------------------
 */
void Direct(double *z_x, double *z_y, double* q, double scale,
            double* double_data, int* int_data, int* realloc_data, 
            int nside, int n_particles, int num_threads) {
    
    /*cursquare is the next box to be treated. Common to all threads.*/
    int i,cursquare=0;
    DirectWorkerStruct *arguments;

    /*The nearest neighbor interaction list.*/
    int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};

    /*Allocate memory for the threads and thread
     *argument structs.*/
    directWorkerThd = malloc(num_threads*sizeof(THREAD_TYPE));
    arguments = malloc(num_threads*sizeof(DirectWorkerStruct));
    
    /*Spawn the threads*/
    for(i=0;i<num_threads;i++) {
        arguments[i].ilist_x      = ilist_x;
        arguments[i].ilist_y      = ilist_y;
        arguments[i].z_x          = z_x;
        arguments[i].z_y          = z_y;
        arguments[i].q            = q;
	arguments[i].scale        = scale;
        arguments[i].double_data  = double_data;
        arguments[i].realloc_data = realloc_data;
        arguments[i].int_data     = int_data;
        arguments[i].n_particles  = n_particles;
        arguments[i].nside        = nside;
        arguments[i].cursquare    = &cursquare;

        THREAD_CREATE(directWorkerThd[i],DirectWorker,(void *) &arguments[i]);
    }

    /*Wait for all threads to complete */
	for(i = 0;i<num_threads;i++)
        THREAD_JOIN(directWorkerThd[i]);

    /*Clean up*/
    free(directWorkerThd);
    free(arguments);

}
/*------------------------------------------------------------------------
 *Threaded worker routine for computing the remaining interactions 
 *directly. Accesses to the output vector is completely parallel, so no
 *mutexes are necessary.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE DirectWorker(void* in_struct) {
    /*Get the structure arrays and variables.*/
    DirectWorkerStruct *arg = (DirectWorkerStruct*)in_struct;
    /*The number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    
    /*Pointer to the output sums.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    /*Pointers to the positions and charges of the particles.*/
    double *z_x = arg->z_x;
    double *z_y = arg->z_y;
    double *q = arg->q;
    double scale = arg->scale;
    double scale2 = scale*scale;
    
    /*The relative offsets of the nearest neighbors.*/
    int *ilist_x = arg->ilist_x;
    int *ilist_y = arg->ilist_y;
    
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    /*The next box to be treated. Common to all threads.*/
    int *cursquare = arg->cursquare;

    /*Loop through the boxes.*/
    for(;;) {
        int i,j,current_box;
        int *tptr;
        
        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&directSqTblmutex);
        i = *cursquare;
        while(i < nside*nside && nparticles_in_box[i] == 0)
            i++;
        if(i==nside*nside) {
            UNLOCK_MUTEX(&directSqTblmutex);
            break;
        }
        current_box = i;
        *cursquare = i+1;
        UNLOCK_MUTEX(&directSqTblmutex);        
        
        /*Temporary pointer to the particles of current box.*/
        tptr = &particle_offsets[box_offsets[current_box]];

        /*Compute the box self-interactions. Since the particles
	 *are in a box, we can compute ij and ji interactions 
	 * simultaneously.*/
        for(j=0;j<nparticles_in_box[current_box];j++) {
            unsigned int k;
            for(k=j+1;k<nparticles_in_box[current_box];k++) {
	      double tz_x = z_x[tptr[k]]-z_x[tptr[j]];
	      double tz_y = z_y[tptr[k]]-z_y[tptr[j]];
	      double tmp = (tz_x*tz_x+tz_y*tz_y);
	      M1_c[tptr[j]] += q[tptr[k]]*expint_log_euler(scale2*tmp);
	      M1_c[tptr[k]] += q[tptr[j]]*expint_log_euler(scale2*tmp);
            }
        }

        /*Compute interactions from the nearest neighbors.*/
        for(j=0;j<8;j++) {
            
            /*Make sure that the box is inside the computational grid.*/
	    if((current_box%nside+ilist_x[j]) < nside && 
               (current_box/nside+ilist_y[j]) < nside &&
                (current_box%nside+ilist_x[j])>= 0 &&
                (current_box/nside+ilist_y[j])>= 0) {

                unsigned int k;
                /*The number of the source nearest neighbor box.*/
                int source_box = current_box+ilist_y[j]*nside+ilist_x[j];

                if(nparticles_in_box[source_box] > 0) {
		  /*Pointer to the particles of the source box.*/
		  int* tptr2 = &particle_offsets[box_offsets[source_box]];
		  
		  /*This loop computes the direct interaction between
		   *all the particles in current_box and the source  
		   *particles. We use SIMD and unroll manually since the 
		   *compiler cannot do that automatically.*/
		  
		  for(k=0;k<(unsigned int)nparticles_in_box[current_box];k++) {
		    unsigned int l=0,l0=0;
		    double res=0;
		    double tz_x, tz_y;
		    double cz_x = z_x[tptr[k]];
		    double cz_y = z_y[tptr[k]];
		    
#ifdef __AVX__
		    double resv[4] __attribute__((aligned(32)));
		    __m256d ZX,ZY,Ql, V;
		    __m256d SCALE2 = _mm256_set1_pd(scale2);
		    __m256d RESV   = _mm256_setzero_pd();
		    __m256d CZX    = _mm256_set1_pd(cz_x);
		    __m256d CZY    = _mm256_set1_pd(cz_y);
	    
		    for(l=0;l<nparticles_in_box[source_box]/8*8;l+=8) {
		      ZX = _mm256_set_pd(z_x[tptr2[l+3]],z_x[tptr2[l+2]],z_x[tptr2[l+1]],z_x[tptr2[l]]);
		      ZY = _mm256_set_pd(z_y[tptr2[l+3]],z_y[tptr2[l+2]],z_y[tptr2[l+1]],z_y[tptr2[l]]);
		      Ql = _mm256_set_pd(q[tptr2[l+3]],q[tptr2[l+2]],q[tptr2[l+1]],q[tptr2[l]]);
		      
		      ZX = _mm256_sub_pd(ZX,CZX);
		      ZX = _mm256_mul_pd(ZX,ZX);
		      ZY = _mm256_sub_pd(ZY,CZY);
		      ZY = _mm256_mul_pd(ZY,ZY);
		      V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
		      V = _mm256_mul_pd(Ql,V);
		      RESV = _mm256_add_pd(RESV, V);
		      
		      l0 = 4;
		      ZX = _mm256_set_pd(z_x[tptr2[l+3+l0]],z_x[tptr2[l+2+l0]],z_x[tptr2[l+1+l0]],z_x[tptr2[l+l0]]);
		      ZY = _mm256_set_pd(z_y[tptr2[l+3+l0]],z_y[tptr2[l+2+l0]],z_y[tptr2[l+1+l0]],z_y[tptr2[l+l0]]);
		      Ql = _mm256_set_pd(q[tptr2[l+3+l0]],q[tptr2[l+2+l0]],q[tptr2[l+1+l0]],q[tptr2[l+l0]]);
		      
		      ZX = _mm256_sub_pd(ZX,CZX);
		      ZX = _mm256_mul_pd(ZX,ZX);
		      ZY = _mm256_sub_pd(ZY,CZY);
		      ZY = _mm256_mul_pd(ZY,ZY);
		      V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
		      V = _mm256_mul_pd(Ql,V);
		      RESV = _mm256_add_pd(RESV, V);
		      }
		    l0 = l;
		    for(l=l0;l<nparticles_in_box[source_box]/4*4;l+=4) {
		      ZX = _mm256_set_pd(z_x[tptr2[l+3]],z_x[tptr2[l+2]],z_x[tptr2[l+1]],z_x[tptr2[l]]);
		      ZY = _mm256_set_pd(z_y[tptr2[l+3]],z_y[tptr2[l+2]],z_y[tptr2[l+1]],z_y[tptr2[l]]);
		      Ql = _mm256_set_pd(q[tptr2[l+3]],q[tptr2[l+2]],q[tptr2[l+1]],q[tptr2[l]]);
		      
		      ZX = _mm256_sub_pd(ZX,CZX);
		      ZX = _mm256_mul_pd(ZX,ZX);
		      ZY = _mm256_sub_pd(ZY,CZY);
		      ZY = _mm256_mul_pd(ZY,ZY);
		      V = _mm256_ein_pd(_mm256_mul_pd(SCALE2,_mm256_add_pd(ZX,ZY)));
		      V = _mm256_mul_pd(Ql,V);
		      RESV = _mm256_add_pd(RESV, V);
		    }
		    //The remainder loop
		    l0 = l;
		    for (l=l0; l<nparticles_in_box[source_box];l++){
		      tz_x = (z_x[tptr2[l]]-cz_x)*(z_x[tptr2[l]]-cz_x);
		      tz_y = (z_y[tptr2[l]]-cz_y)*(z_y[tptr2[l]]-cz_y);
		      res += q[tptr2[l]]*expint_log_euler( scale2*(tz_x+tz_y) );
		    }
		  
		    _mm256_store_pd(resv,RESV);
		    res += sum4(resv);

#elif defined __SSE4_2__
		    double resv[2] MEM_ALIGN;
		    __m128d ZX,ZY,Ql, V;
		    __m128d SCALE2 = _mm_set1_pd(scale2);
		    __m128d RESV = _mm_setzero_pd();
		    __m128d CZX = _mm_set1_pd(cz_x);
		    __m128d CZY = _mm_set1_pd(cz_y);
		    for(l=0;l<nparticles_in_box[source_box]/4*4;l+=4) {
		      ZX = _mm_set_pd(z_x[tptr2[l+1]],z_x[tptr2[l]]);
		      ZY = _mm_set_pd(z_y[tptr2[l+1]],z_y[tptr2[l]]);
		      Ql = _mm_set_pd(q[tptr2[l+1]],q[tptr2[l]]);
		      
		      ZX = _mm_sub_pd(ZX,CZX);
		      ZX = _mm_mul_pd(ZX,ZX);
		      ZY = _mm_sub_pd(ZY,CZY);
		      ZY = _mm_mul_pd(ZY,ZY);
		      V = _mm_ein_pd(_mm_mul_pd(SCALE2,_mm_add_pd(ZX,ZY)));
		      V = _mm_mul_pd(Ql,V);
		      RESV = _mm_add_pd(RESV, V);
		      
		      l0 = 2;
		      ZX = _mm_set_pd(z_x[tptr2[l+1+l0]],z_x[tptr2[l+l0]]);
		      ZY = _mm_set_pd(z_y[tptr2[l+1+l0]],z_y[tptr2[l+l0]]);
		      Ql = _mm_set_pd(q[tptr2[l+1+l0]],q[tptr2[l+l0]]);
		      
		      ZX = _mm_sub_pd(ZX,CZX);
		      ZX = _mm_mul_pd(ZX,ZX);
		      ZY = _mm_sub_pd(ZY,CZY);
		      ZY = _mm_mul_pd(ZY,ZY);
		      V = _mm_ein_pd(_mm_mul_pd(SCALE2,_mm_add_pd(ZX,ZY)));
		      V = _mm_mul_pd(Ql,V);
		      RESV = _mm_add_pd(RESV, V);
		    }
		    //The remainder loop
		    l0 = l;
		    for (l=l0; l<nparticles_in_box[source_box];l++){
		      tz_x = (z_x[tptr2[l]]-cz_x)*(z_x[tptr2[l]]-cz_x);
		      tz_y = (z_y[tptr2[l]]-cz_y)*(z_y[tptr2[l]]-cz_y);
		      res += q[tptr2[l]]*expint_log_euler( scale2*(tz_x+tz_y));
		    }
		    
		    _mm_store_pd(resv,RESV);
		    res += resv[0]+resv[1];
#else

		    for (l=0; l<nparticles_in_box[source_box];l++){
		      tz_x = (z_x[tptr2[l]]-cz_x)*(z_x[tptr2[l]]-cz_x);
		      tz_y = (z_y[tptr2[l]]-cz_y)*(z_y[tptr2[l]]-cz_y);
		      res += q[tptr2[l]]*expint_log_euler( scale2*(tz_x+tz_y) );
		    }
#endif
		    M1_c[tptr[k]] += res;
		  }
		}
	    }
	}
    }
    THREAD_EXIT();
}

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *z_x, double *z_y, int n_particles, int nside,
            int *maxparticles_in_box, int *int_data, int *realloc_data) {
            
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int* particle_offsets = &int_data[PARTICLE_OFFSETS_OFFSET];
    /*The starting offsets in particle_offsets for each box.*/
    int* box_offsets = &realloc_data[BOX_OFFSETS_OFFSET];
    /*The number of particles in each box.*/
    int* nparticles_in_box = &realloc_data[NPARTICLES_IN_BOX_OFFSET];
    /*Temporary arrays.*/
    int* in_box = &int_data[IN_BOX_OFFSET];
    int* particlenum_in_box = &realloc_data[ASSIGN_TMP_OFFSET];
    /*The total number of boxes.*/
    int number_of_boxes = nside*nside;
    int i,j;
    
    /*Clear the nparticles_in_box array. It is realloc:d for each grid.*/
    for(j = 0;j<number_of_boxes;j++){
        nparticles_in_box[j] = 0;
    }

    /*Assign the particles to boxes.*/
    for(j = 0;j<n_particles;j++) {
        int box_x = (int)floor(nside*(z_x[j]+.5));
        int box_y = (int)floor(nside*(z_y[j]+.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside) box_x = nside-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside) box_y = nside-1;
        in_box[j] =  box_y*nside + box_x;
        nparticles_in_box[in_box[j]]++;
    }

    /*Keep track of the maximum number of particles in a box. Wen this
     *quantity falls below n_terms, we do one last sweep of multipole
     *calculations and then compute the rest of the interactions directly.
     */
    maxparticles_in_box[0] = -1;
    for(j=0;j<number_of_boxes;j++) {
        if(nparticles_in_box[j] > *maxparticles_in_box)
            *maxparticles_in_box = nparticles_in_box[j];
    }
    
    /*box_offsets is the offsets of the boxes in particle_offsets to be 
     *computed below. particlenum_in_box is a temporary array.*/
    box_offsets[0] = particlenum_in_box[0] = 0;
    i = 0;
    for(j=1;j<number_of_boxes+1;j++) {
        i += nparticles_in_box[j-1];
        particlenum_in_box[j] = box_offsets[j] = i;
    }
    
    /*The offsets of the particles in each box in the z and q arrays.*/
    for(j=0;j<n_particles;j++) 
        particle_offsets[particlenum_in_box[in_box[j]]++] = j;
}


/*------------------------------------------------------------------------
 *This function computes the multipole coefficients.
 *------------------------------------------------------------------------
 */
void compute_mpole_c(double* mpole_c, int nt, double t_x, double t_y)
{
  int k1, k2, n;
  double x2y2 = t_x*t_x + t_y*t_y;
  double exp2 = exp(-x2y2);
  double s=0;

  double* mpole_b = _mm_mxMalloc(nt*nt*sizeof(double),32);

  /* compute mpole_b coeffs*/
  mpole_b[0*nt+0] =  exp2;
  mpole_b[1*nt+0] = -2.0*t_x*exp2;
  mpole_b[0*nt+1] = -2.0*t_y*exp2;
  mpole_b[1*nt+1] =  4.0*t_x*t_y*exp2;

  for (k1=2; k1<nt; k1++){
    double Neg_Two_k1 = -2.0/(double) k1;
    mpole_b[k1*nt+0] = Neg_Two_k1 *( t_x*mpole_b[(k1-1)*nt+0] + mpole_b[(k1-2)*nt+0] );
    mpole_b[k1*nt+1] = Neg_Two_k1 *( t_x*mpole_b[(k1-1)*nt+1] + mpole_b[(k1-2)*nt+1] );
    mpole_b[0*nt+k1] = Neg_Two_k1 *( t_y*mpole_b[0*nt+(k1-1)] + mpole_b[0*nt+(k1-2)] );
    mpole_b[1*nt+k1] = Neg_Two_k1 *( t_y*mpole_b[1*nt+(k1-1)] + mpole_b[1*nt+(k1-2)] );
  }

  if(nt>=4){
    /* for (k1=2; k1<nt; k1++) */
    /*   for (k2=2; k2<nt; k2++){ */
    for (n=4; n<nt; n++)
      for(k1=2;k1<=(n-2); k1++){
	k2 = n-k1;
	double Neg_Two_k2 = -2./(double) k2;
	mpole_b[k1*nt+k2] = Neg_Two_k2 * ( t_y*mpole_b[k1*nt+(k2-1)] + mpole_b[k1*nt+(k2-2)] );
      }
  }
    /* This accounts for the presence of the log term */
  mpole_b[0*nt+0] += x2y2;
  mpole_b[1*nt+0] += 2.0*t_x;
  mpole_b[0*nt+1] += 2.0*t_y;
  mpole_b[1*nt+1] += 0.0;
  mpole_b[2*nt+0] += 1.0;
  mpole_b[0*nt+2] += 1.0;
  
  /* compute mpole_a coeffs*/
  mpole_c[0*nt+0] = expint_log_euler(x2y2);
  mpole_c[1*nt+0] = 2.0*t_x*(1.0-exp2)/x2y2;
  mpole_c[0*nt+1] = 2.0*t_y*(1.0-exp2)/x2y2;
  mpole_c[1*nt+1] = 4.0*t_x*t_y*(exp2*(x2y2+1.0) - 1.0)/(x2y2*x2y2);
                            
  /*compute the rest of mpole_a coeffs and add up to form mpole_c*/
  for (k1=2; k1<nt; k1++){
    k2 = 0; s = (double) k1+k2;
    mpole_c[k1*nt+k2] = (-2.*(s-1.0)*t_x*mpole_c[(k1-1)*nt+0]
			 -(s-2.0)*mpole_c[(k1-2)*nt+0] + mpole_b[k1*nt+0]*s)/s/x2y2;
    mpole_c[k2*nt+k1] = (-2.*(s-1.0)*t_y*mpole_c[0*nt+(k1-1)]
			 -(s-2.0)*mpole_c[0*nt+(k1-2)] + mpole_b[0*nt+k1]*s)/s/x2y2;
    k2 = 1; s = (double) k1+k2;
    mpole_c[k1*nt+k2] = (-2.*(s-1.)* (t_x*mpole_c[(k1-1)*nt+k2] + t_y*mpole_c[k1*nt+(k2-1)]) 
			 -(s-2.)*mpole_c[(k1-2)*nt+k2] + mpole_b[k1*nt+k2]*s)/s/x2y2;
    mpole_c[k2*nt+k1] = (-2.*(s-1.)* (t_x*mpole_c[(k2-1)*nt+k1] + t_y*mpole_c[k2*nt+(k1-1)]) 
			 -(s-2.)*mpole_c[k2*nt+(k1-2)] + mpole_b[k2*nt+k1]*s)/s/x2y2;
  }


  if(nt>=4){
  /* for (k1=2; k1<nt; k1++) */
  /*   for (k2=2; k2<nt; k2++){ */
    for (n=4; n<nt; n++){
      for(k1=2;k1<=(n-2); k1++){
	k2 = n-k1;
	s = (double) k1+k2;
	mpole_c[k1*nt+k2] = (-2.*(s-1.) * (t_x*mpole_c[(k1-1)*nt+k2]+t_y*mpole_c[k1*nt+(k2-1)])
			     -(s-2.) * (    mpole_c[(k1-2)*nt+k2]+    mpole_c[k1*nt+(k2-2)])
			     + s     *      mpole_b[k1*nt+k2]   )/s/x2y2;
      }
    }
  }
  _mm_mxFree(mpole_b);
}
