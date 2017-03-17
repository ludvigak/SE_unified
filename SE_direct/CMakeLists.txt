## FGG MEX libraries
  matlab_add_mex(
    NAME SE3P_direct_real_mex
    SRC SE3P_direct_real.c
    )
  matlab_add_mex(
    NAME SE3P_direct_rsrc_mex
    SRC SE3P_direct_rsrc.c
    ) 
  matlab_add_mex(
    NAME SE3P_direct_fd_mex
    SRC SE3P_direct_fd.c
    )  
  matlab_add_mex(
    NAME SE3P_direct_self_mex
    SRC SE3P_direct_self.c
    )
  
  matlab_add_mex(
    NAME SE2P_direct_real_mex
    SRC SE2P_direct_real.c
    ) 
  matlab_add_mex(
    NAME SE2P_direct_rsrc_mex
    SRC SE2P_direct_rsrc.c
    ) 
  matlab_add_mex(
    NAME SE2P_direct_fd_mex
    SRC SE2P_direct_fd.c
    )
  matlab_add_mex(
    NAME SE2P_direct_self_mex
    SRC SE2P_direct_self.c
    )
  matlab_add_mex(
    NAME SE2P_direct_k0_mex
    SRC SE2P_direct_k0.c
    )

