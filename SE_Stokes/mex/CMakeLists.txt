
add_definitions(-DHASIMOTO -DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

matlab_add_mex(
  NAME stokeslet_fast_k_scaling
  SRC stokeslet_fast_k_scaling.c
)
