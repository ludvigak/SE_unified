#include_directories(
#  ${CMAKE_CURRENT_SOURCE_DIR}
#)

add_definitions(-DBEENAKKER -DMATLAB_MEX_FILE)

matlab_add_mex(
  NAME stresslet_fast_k_scaling
  SRC stresslet_fast_k_scaling.c
)

matlab_add_mex(
  NAME stresslet_direct_real_mexcore
  SRC stresslet_direct_real_mexcore.c
)

matlab_add_mex(
  NAME stresslet_real_rc_mex
  SRC stresslet_real_rc_mex.c stresslet_real_rc.c cell_list.c
)

matlab_add_mex(
  NAME stresslet_real_rc_nomatrix_mex
  SRC stresslet_real_rc.c cell_list.c stresslet_real_rc_nomatrix_mex.c
)
# largeArrayDims -DBEENAKKER  mex/stresslet_real_rc_mex