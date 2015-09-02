#ifndef __FGG_THRD_H
#define __FGG_THRD_H

#include "omp.h"

typedef struct {
    int p;    
    int block_begin;
    int block_end;
    int block_size;
    int skip;
} grid_thrd_ws_t;

// If FGG_THRD activated without OpenMP, throw warning and disable
#ifdef FGG_THRD
#ifndef _OPENMP
#warning "FGG_THRD will not work without OpenMP"
#undef FGG_THRD
#endif
#endif

#endif
