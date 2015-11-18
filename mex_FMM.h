/*
 * Header file containing definitions for mex:ed fast multipole
 */

#ifndef _MEX_FMM_H
#define _MEX_FMM_H

#include <math.h>
#include <sys/time.h>

/*Required to utilize SSE-intrinsics*/
#include <immintrin.h>

/*Contains functions to allocate aligned memory using Matlabs
 *mxMalloc and mxCalloc*/
#include "mm_mxmalloc.h"

/*We use threads. Try to make this at least somewhat portable, so that
 *we can compile and run on LINUX/MAC OS X as well as on Windows.*/
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

/*Only needed to compute the required number of multipole terms*/
#define LOG2 0.69314718055994529

/*The size of the preallocation of some of the data, before we resort to
 *reallocation. This needs to be an even power of two.*/
#define PREALLOCNSIDE (128)
/*Flags that a box is treated by a thread*/
#define TREATED (2100000000)

/*Some error messages*/
#define Z_CLASS_ERROR "z must be real or complex."
#define Q_CLASS_ERROR "q must be real or complex."
#define TOL_CLASS_ERROR "tol must be a real scalar."
#define LEVMAX_CLASS_ERROR "levmax must be a real scalar."
#define LEVMAX_WARNING "levmax not allowed to be >16. Setting levmax = 16."
#define NUMTHREADS_CLASS_ERROR "numthreads must be a real scalar."
#define Z_DIM_ERROR "z must be a column-vector."
#define Q_DIM_ERROR "q must be a column-vector."
#define Z_Q_SIZE_ERROR "length(z) is not equal to length(q)."
#define OUTSIDE_ERROR "Points outside computational box([-0.5,0.5,-0.5,0.5])."
#define NUM_THREADS_TOO_LOW_ERROR "At least one thread required."

/*Define thread related macros and types. This is a crude wrapper to make
 *the code somewhat portable.*/
#ifdef WIN32
typedef CRITICAL_SECTION MUTEX_TYPE; 
typedef HANDLE THREAD_TYPE; 
/*Ugly!*/
#define THREAD_FUNC_TYPE DWORD WINAPI 
#define INIT_MUTEX(MUTEX) InitializeCriticalSection(MUTEX)
#define DESTROY_MUTEX(MUTEX) DeleteCriticalSection(MUTEX)
#define LOCK_MUTEX(MUTEX) EnterCriticalSection(MUTEX)
#define UNLOCK_MUTEX(MUTEX) LeaveCriticalSection(MUTEX)
#define THREAD_CREATE(THREAD,FUNC,ARGS) THREAD=CreateThread(NULL, 0, FUNC, ARGS, 0, NULL)
#define THREAD_JOIN(THREAD) WaitForSingleObject(THREAD,INFINITE)
#define THREAD_EXIT() return 0
#else
typedef pthread_mutex_t MUTEX_TYPE;
typedef pthread_t THREAD_TYPE; 
/*Ugly!*/
#define THREAD_FUNC_TYPE void*
#define INIT_MUTEX(MUTEX) pthread_mutex_init(MUTEX,NULL)
#define DESTROY_MUTEX(MUTEX) pthread_mutex_destroy(MUTEX)
#define LOCK_MUTEX(MUTEX) pthread_mutex_lock(MUTEX)
#define UNLOCK_MUTEX(MUTEX) pthread_mutex_unlock(MUTEX)
#define THREAD_CREATE(THREAD,FUNC,ARGS) pthread_create(&THREAD, NULL, FUNC, ARGS)
#define THREAD_JOIN(THREAD) pthread_join(THREAD,NULL)
#define THREAD_EXIT() pthread_exit((void*) 0)
#endif

/*The threads and mutexes. The number of mutexes depends on the number of 
 *threads wanted.*/
int num_mutexes;
THREAD_TYPE *mpoleWorkerThd,*directWorkerThd;
MUTEX_TYPE   mpoleSqTblmutex,directSqTblmutex;
MUTEX_TYPE  *output_mutex,*localexp_mutex;

/*The thread worker functions.*/
/*------------------------------------------------------------------------
 *Threaded worker function that computes multipole expansions for the
 *particles in each box and then traverses the interaction list doing one
 *of three things: either the interactions are computed directly, or the 
 *multipole expansion is evaluated directly in the target box, or the 
 *multipole expansions are converted into local taylor expansions around 
 *the centers of the target boxes and summed up.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorker(void* in_struct);
/*------------------------------------------------------------------------
 *Threaded worker function that evaluates and sums the local taylor 
 *expansions that were computed in MpolesWorker.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorkerSum(void* in_struct);
/*------------------------------------------------------------------------
 *Threaded worker routine for computing the remaining interactions 
 *directly.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE DirectWorker(void* in_struct);


/*------------------------------------------------------------------------
 * Initialize the interaction list. Same for all multipole codes.
 *------------------------------------------------------------------------
 */
void Intlist(int* ilist_x,int* ilist_y) {
    int i;
    for(i=0;i<4;i++) {
        ilist_x[0+i*27]=ilist_x[4+i*27]=ilist_x[8+i*27]=ilist_x[11+i*27]=ilist_x[15+i*27]=2;
        ilist_x[2+i*27]=ilist_x[6+i*27]=ilist_x[9+i*27]=ilist_x[10+i*27]=ilist_x[13+i*27]=-2;
        ilist_x[5+i*27]=ilist_x[14+i*27]=-1;
        ilist_x[7+i*27]=ilist_x[12+i*27]=1;
        ilist_y[1+i*27]=ilist_y[5+i*27]=ilist_y[8+i*27]=ilist_y[9+i*27]=ilist_y[12+i*27]=2;
        ilist_y[3+i*27]=ilist_y[7+i*27]=ilist_y[10+i*27]=ilist_y[11+i*27]=ilist_y[14+i*27]=-2;
        ilist_y[4+i*27]=ilist_y[13+i*27]=1;
        ilist_y[6+i*27]=ilist_y[15+i*27]=-1;
    }
    ilist_x[16]=ilist_x[18]=ilist_x[20]=ilist_x[22]=ilist_x[24]=ilist_x[26]=
            ilist_x[44]=ilist_x[46]=ilist_x[48]=ilist_x[50]=ilist_x[52]=ilist_x[53]=3;
    ilist_x[17]=ilist_x[51]=ilist_x[70]=ilist_x[106]=2;
    ilist_x[19]=ilist_x[49]=ilist_x[72]=ilist_x[104]=1;
    ilist_x[23]=ilist_x[45]=ilist_x[76]=ilist_x[100]=-1;
    ilist_x[25]=ilist_x[43]=ilist_x[78]=ilist_x[98]=-2;
    ilist_x[71]=ilist_x[73]=ilist_x[75]=ilist_x[77]=ilist_x[79]=ilist_x[80]=
            ilist_x[97]=ilist_x[99]=ilist_x[101]=ilist_x[103]=ilist_x[105]=ilist_x[107]=-3;
    
    ilist_y[17]=ilist_y[19]=ilist_y[21]=ilist_y[23]=ilist_y[25]=ilist_y[26]=
            ilist_y[70]=ilist_y[72]=ilist_y[74]=ilist_y[76]=ilist_y[78]=ilist_y[80]=3;
    ilist_y[24]=ilist_y[52]=ilist_y[71]=ilist_y[97]=2;
    ilist_y[22]=ilist_y[50]=ilist_y[73]=ilist_y[99]=1;
    ilist_y[18]=ilist_y[46]=ilist_y[77]=ilist_y[103]=-1;
    ilist_y[16]=ilist_y[44]=ilist_y[79]=ilist_y[105]=-2;
    ilist_y[43]=ilist_y[45]=ilist_y[47]=ilist_y[49]=ilist_y[51]=ilist_y[53]=
            ilist_y[98]=ilist_y[100]=ilist_y[102]=ilist_y[104]=ilist_y[106]=ilist_y[107]=-3;
}


/*------------------------------------------------------------------------
 * Time the code
 *------------------------------------------------------------------------
 */
double gettime(void)
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

void compute_mpole_c(double*, int, double, double);

void isaligned(void* p, unsigned int n)
{
    if ( ((unsigned long) p & (n-1))==0)
        printf("is %d bit aligned\n", n);
    else
        printf("is NOT %d bit aligned\n",n);
}

unsigned int MPOLE=0, TAYLOR=0, DIRECT=0;

typedef struct T{
  double START;
  double STOP;
} TIME0;

TIME0 TIME;

void START_TIME()
{
   TIME.START = gettime();
}
void STOP_TIME()
{
   TIME.STOP += gettime()-TIME.START;
}

void print_log()
{
  printf("RUNTIME %f\n",TIME.STOP);
}


/* sum up 4 elements in a vector, used for SIMD op.*/
double sum4(double *p)
{
  return (p[0]+p[1]+p[2]+p[3]);
}

#endif
