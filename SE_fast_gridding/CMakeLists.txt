include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

## FGG MEX libraries
add_library(
  fgg
  SE_fgg.c SE_fgg.h SE_fgg_MEX_params.c
)
target_compile_definitions(fgg PUBLIC THREE_PERIODIC)

add_library(
  fgg_thrd
  SE_fgg.c SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd PUBLIC THREE_PERIODIC FGG_THRD)

## MEX functions
matlab_add_mex(
  NAME SE_fg_grid_mex
  SRC mex/SE_fg_grid_mex.c 
)
target_link_libraries(SE_fg_grid_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_mex
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex fgg)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex
  SRC mex/SE_fgg_expand_all_mex.c  SE_fgg.c SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex fgg)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex fgg)

matlab_add_mex(
  NAME SE_fg_grid_split_mex
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_split_mex
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex fgg)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex fgg_thrd)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex fgg)
