## FGG MEX libraries
add_library(
  mathint
  mathint.c mathint.h SE_direct.h math_x86.h
)
target_compile_definitions(mathint PUBLIC POTENTIAL gls gslcblas)

matlab_add_mex(
  NAME SE1P_direct_mex
  SRC SE1P_direct.c
  )
target_link_libraries(SE1P_direct_mex mathint)

matlab_add_mex(
  NAME SE1P_direct_real_mex
  SRC SE1P_direct_real.c
  )
target_link_libraries(SE1P_direct_real_mex mathint)

matlab_add_mex(
  NAME SE1P_direct_rsrc_mex
  SRC SE1P_direct_rsrc.c
  )
target_link_libraries(SE1P_direct_rsrc_mex mathint)

matlab_add_mex(
  NAME SE1P_direct_fd_mex
  SRC SE1P_direct_fd.c # mathint.h
  )
target_link_libraries(SE1P_direct_fd_mex mathint gsl gslcblas)

matlab_add_mex(
  NAME SE1P_direct_k0_mex
  SRC SE1P_direct_k0.c
  )
target_link_libraries(SE1P_direct_k0_mex mathint gsl gslcblas)

matlab_add_mex(
  NAME SE1P_direct_self_mex
  SRC SE1P_direct_self.c
  )
target_link_libraries(SE1P_direct_self_mex mathint)


add_library(
  mathint_force
  mathint.c mathint.h SE_direct.h math_x86.h
)
target_compile_definitions(mathint_force PUBLIC FORCE)

matlab_add_mex(
  NAME SE1P_direct_force_mex
  SRC SE1P_direct.c
  )
target_link_libraries(SE1P_direct_force_mex mathint_force)

matlab_add_mex(
  NAME SE1P_direct_real_force_mex
  SRC SE1P_direct_real.c
  )
target_link_libraries(SE1P_direct_real_force_mex mathint_force)

matlab_add_mex(
  NAME SE1P_direct_rsrc_force_mex
  SRC SE1P_direct_rsrc.c
  )
target_link_libraries(SE1P_direct_rsrc_force_mex mathint_force)

matlab_add_mex(
  NAME SE1P_direct_fd_force_mex
  SRC SE1P_direct_fd.c
  )
target_link_libraries(SE1P_direct_fd_force_mex mathint_force gsl gslcblas)

matlab_add_mex(
  NAME SE1P_direct_k0_force_mex
  SRC SE1P_direct_k0.c
  )
target_link_libraries(SE1P_direct_k0_force_mex mathint_force)
