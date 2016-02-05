
add_definitions(-DHASIMOTO -DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

## 2P
# 2P real
matlab_add_mex(
  NAME SE2P_Stokes_direct_real_mex
  SRC SE2P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE2P_Stokes_direct_real_mex PUBLIC TWO_PERIODIC EWALD_REAL)

# 2P rsrc
matlab_add_mex(
  NAME SE2P_Stokes_direct_rsrc_mex
  SRC SE2P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE2P_Stokes_direct_rsrc_mex PUBLIC TWO_PERIODIC EWALD_RSRC)

# 2P fd
matlab_add_mex(
  NAME SE2P_Stokes_direct_fd_mex
  SRC SE2P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE2P_Stokes_direct_fd_mex PUBLIC TWO_PERIODIC EWALD_FD)

# 2P k0
matlab_add_mex(
  NAME SE2P_Stokes_direct_k0_mex
  SRC SE2P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE2P_Stokes_direct_k0_mex PUBLIC TWO_PERIODIC EWALD_K0)

# 2P self
matlab_add_mex(
  NAME SE2P_Stokes_direct_self_mex
  SRC SE2P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE2P_Stokes_direct_self_mex PUBLIC TWO_PERIODIC EWALD_SELF)

## 3P
# 3P real
matlab_add_mex(
  NAME SE3P_Stokes_direct_real_mex
  SRC SE3P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE3P_Stokes_direct_real_mex PUBLIC THREE_PERIODIC EWALD_REAL)

# 3P rsrc
matlab_add_mex(
  NAME SE3P_Stokes_direct_rsrc_mex
  SRC SE3P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE3P_Stokes_direct_rsrc_mex PUBLIC THREE_PERIODIC EWALD_RSRC)

# 3P rsrc ext
matlab_add_mex(
  NAME SE3P_Stokes_direct_rsrc_ext_mex
  SRC SE3P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE3P_Stokes_direct_rsrc_ext_mex PUBLIC THREE_PERIODIC EWALD_RSRC EXTERNAL)


# 3P fd
matlab_add_mex(
  NAME SE3P_Stokes_direct_fd_mex
  SRC SE3P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE3P_Stokes_direct_fd_mex PUBLIC THREE_PERIODIC EWALD_FD)

# 3P self
matlab_add_mex(
  NAME SE3P_Stokes_direct_self_mex
  SRC SE3P_Stokes_direct.c SE_Stokes_direct_mex.c
)
target_compile_definitions(SE3P_Stokes_direct_self_mex PUBLIC THREE_PERIODIC EWALD_SELF)




