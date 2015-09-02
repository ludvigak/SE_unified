#ifndef __SE_FGG_THRD_H
#define __SE_FGG_THRD_H

#include "SE_fgg.h"

void SE_FGG_grid_thrd(SE_FGG_work* work, const SE_state* st, 
		     const SE_FGG_params* params);

// If FGG_THRD activated without OpenMP, throw warning and disable
#ifdef FGG_THRD
#ifndef _OPENMP
#warning "FGG_THRD will not work without OpenMP"
#undef FGG_THRD
#endif
#endif

#endif
