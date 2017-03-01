#add_definitions(-DVERBOSE )
add_definitions(-DHASIMOTO )

# 2P kspace
matlab_add_mex(
  NAME mul_k_space_kernel_mex
  SRC mul_k_space_kernel_mex.c
)
