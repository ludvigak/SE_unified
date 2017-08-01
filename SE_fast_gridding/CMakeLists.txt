include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

## FGG MEX libraries
# 3D
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

add_library(
  fgg_force
  SE_fgg.c SE_fgg.h fgg_force.c SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_force PUBLIC THREE_PERIODIC FORCE)

## FGG MEX libraries

# 2D
add_library(
  fgg_2p
  SE_fgg.c SE_fgg.h SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_2p PUBLIC TWO_PERIODIC)

add_library(
  fgg_thrd_2p
  SE_fgg.c SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd_2p PUBLIC TWO_PERIODIC FGG_THRD)

## FGG MEX libraries

# 1D

add_library(
  fgg_1p
  SE_fgg.c SE_fgg.h SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_1p PUBLIC ONE_PERIODIC)

add_library(
  fgg_thrd_1p
  SE_fgg.c SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd_1p PUBLIC ONE_PERIODIC FGG_THRD)

add_library(
  fgg_force_1p
  SE_fgg.c SE_fgg.h SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_force_1p PUBLIC ONE_PERIODIC FORCE)


## MEX functions 3P

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

## MEX functions 2P
matlab_add_mex(
  NAME SE_fg_grid_mex_2p
  SRC mex/SE_fg_grid_mex.c 
)
target_link_libraries(SE_fg_grid_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_mex_2p
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex_2p
  SRC mex/SE_fgg_expand_all_mex.c  SE_fgg.c SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex_2p
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_mex_2p
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_split_mex_2p
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex_2p
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex_2p fgg_thrd_2p)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex_2p
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex_2p
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex_2p fgg_2p)

## MEX functions 1P
matlab_add_mex(
  NAME SE_fg_grid_mex_1p
  SRC mex/SE_fg_grid_mex.c 
)
target_link_libraries(SE_fg_grid_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_mex_1p
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex_1p
  SRC mex/SE_fgg_expand_all_mex.c  SE_fgg.c SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex_1p
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_mex_1p
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_split_mex_1p
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex_1p
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex_1p fgg_thrd_1p)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex_1p
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex_1p
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex_1p fgg_1p)

# MEX functions 1P FORCE
matlab_add_mex(
  NAME SE_fgg_expand_all_force_mex_1p
  SRC mex/SE_fgg_expand_all_force_mex.c  SE_fgg.c SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_force_mex_1p fgg_force_1p)

matlab_add_mex(
  NAME SE_fg_int_force_mex_1p
  SRC mex/SE_fg_int_force_mex.c
)
target_link_libraries(SE_fg_int_force_mex_1p fgg_force_1p)

matlab_add_mex(
  NAME SE_fg_int_split_force_mex_1p
  SRC mex/SE_fg_int_split_force_mex.c
)
target_link_libraries(SE_fg_int_split_force_mex_1p fgg_force_1p)


############## KAISER MEX FUNCTIONS
## FGG MEX libraries
# 3D
add_library(
  fkg
  SE_fgg.c SE_fgg.h SE_fkg_MEX_params.c
)
target_compile_definitions(fkg PUBLIC THREE_PERIODIC)

## FGG MEX libraries

# 2D
add_library(
  fkg_2p
  SE_fgg.c SE_fgg.h SE_fkg_MEX_params.c
)
target_compile_definitions(fkg_2p PUBLIC TWO_PERIODIC)

## FGG MEX libraries

# 1D

add_library(
  fkg_1p
  SE_fgg.c SE_fgg.h SE_fkg_MEX_params.c
)
target_compile_definitions(fkg_1p PUBLIC ONE_PERIODIC)


#  KAISER MEX functions 3P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex
  SRC mex/SE_fg_grid_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex
  SRC mex/SE_fg_grid_split_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_split_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fg_extend_fcn_kaiser_mex
  SRC mex/SE_fg_extend_fcn_kaiser_mex.c
)
target_link_libraries(SE_fg_extend_fcn_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex fkg)

#  KAISER MEX functions 2P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex_2p
  SRC mex/SE_fg_grid_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_kaiser_mex_2p fkg_2p)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex_2p
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex_2p fkg_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex_2p
  SRC mex/SE_fg_grid_split_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_split_kaiser_mex_2p fkg_2p)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex_2p
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex_2p fkg_2p)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex_2p
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex_2p fkg_2p)

#  KAISER MEX functions 1P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex_1p
  SRC mex/SE_fg_grid_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_kaiser_mex_1p fkg_1p)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex_1p
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex_1p fkg_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex_1p
  SRC mex/SE_fg_grid_split_kaiser_mex.c 
)
target_link_libraries(SE_fg_grid_split_kaiser_mex_1p fkg_1p)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex_1p
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex_1p fkg_1p)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex_1p
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex_1p fkg_1p)
