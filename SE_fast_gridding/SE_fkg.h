#ifndef __SE_FKG_H
#define __SE_FKG_H

#include "SE_fgg.h"

// kaiser signatures
void SE_FKG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int);
void SE_FKG_expand_all(SE_FGG_work*, const SE_state* , const SE_FGG_params*);
void SE_FKG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int(double*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
/* void SE_FKG_grid_split(SE_FGG_work*, const SE_state*, const SE_FGG_params*); */
/* void SE_FKG_int_split(double*, const SE_FGG_work*, const SE_FGG_params*); */
/* void SE_FKG_int_split_SSE(double*, const SE_FGG_work*, const SE_FGG_params*); */
/* void SE_FKG_int_split_SSE_P8(double*, const SE_FGG_work*, const SE_FGG_params*); */
/* void SE_FKG_int_split_AVX_P8(double*, const SE_FGG_work*, const SE_FGG_params*); */
/* void SE_FKG_grid_split_SSE(SE_FGG_work*, const SE_state*, const SE_FGG_params*); */
/* void SE_FKG_grid_split_SSE_u8(SE_FGG_work*, const SE_state*, const SE_FGG_params*); */
/* void SE_FKG_grid_split_AVX_P8(SE_FGG_work*, const SE_state*, const SE_FGG_params*); */

void SE_FKG_int_split_SSE_dispatch(double*, const SE_FGG_work*, const SE_FGG_params* );
void SE_FKG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int_split_AVX_dispatch(double*, const SE_FGG_work*, const SE_FGG_params* );
void SE_FKG_grid_split_AVX_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
#endif
