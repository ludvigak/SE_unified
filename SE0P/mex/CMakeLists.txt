## MEX functions

matlab_add_mex(
  NAME SE0P_stokeslet_fast_fs_k_scaling
  SRC stokeslet_fast_fs_k_scaling.c
)

matlab_add_mex(
  NAME SE0P_stokeslet_direct_mex
  SRC stokeslet_direct.c
)

matlab_add_mex(
  NAME SE0P_stokeslet_rsrc_mex
  SRC stokeslet_rsrc.c
)

matlab_add_mex(
  NAME SE0P_stresslet_fast_fs_k_scaling
  SRC stresslet_fast_fs_k_scaling.c
)

matlab_add_mex(
  NAME SE0P_stresslet_direct_mex
  SRC stresslet_direct.c
)

matlab_add_mex(
  NAME SE0P_stresslet_rsrc_mex
  SRC stresslet_rsrc.c
)

matlab_add_mex(
  NAME SE0P_rotlet_direct_mex
  SRC rotlet_direct.c
)

matlab_add_mex(
  NAME SE0P_rotlet_rsrc_mex
  SRC rotlet_rsrc.c
)

if(MKL_FOUND)
  
  # Build stokeslet RS using cell list and MKL
  matlab_add_mex(
    NAME SE0P_stokeslet_rsrc_cell_mkl_mex
    SRC stokeslet_rsrc_cell.c cell_list.c
    )
  
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE0P_stokeslet_rsrc_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE0P_stokeslet_rsrc_cell_mkl_mex ${MKL_STATIC})

  target_compile_definitions(SE0P_stokeslet_rsrc_cell_mkl_mex PUBLIC INTEL_MKL)
  
endif()

# no-MKL version
matlab_add_mex(
  NAME SE0P_stokeslet_rsrc_cell_mex
  SRC stokeslet_rsrc_cell.c cell_list.c
)
