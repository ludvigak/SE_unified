include(FindPackageHandleStandardArgs)

find_path(MKL_ROOT include/mkl.h
  PATHS
  /opt/intel/mkl
  $ENV{MKLROOT}
  $ENV{HOME}/intel/mkl
)

find_path(MKL_INCLUDE_DIR mkl.h
  PATHS
  ${MKL_ROOT}/include  
)

find_path(MKL_LIBRARY libmkl_core.a
  PATHS
  ${MKL_ROOT}/lib/intel64
)

# List of all libraries
set(MKL_ALL_STATIC
  libmkl_blacs_intelmpi_ilp64.a  libmkl_cdft_core.a     libmkl_lapack95_ilp64.a
  libmkl_blacs_intelmpi_lp64.a   libmkl_core.a          libmkl_lapack95_lp64.a
  libmkl_blacs_openmpi_ilp64.a   libmkl_gf_ilp64.a      libmkl_scalapack_ilp64.a
  libmkl_blacs_openmpi_lp64.a    libmkl_gf_lp64.a       libmkl_scalapack_lp64.a
  libmkl_blacs_sgimpt_ilp64.a    libmkl_gnu_thread.a    libmkl_sequential.a
  libmkl_blacs_sgimpt_lp64.a     libmkl_intel_ilp64.a   libmkl_tbb_thread.a
  libmkl_blas95_ilp64.a          libmkl_intel_lp64.a
  libmkl_blas95_lp64.a           libmkl_intel_thread.a
)

# Setup intel MKL
set(MKL_STATIC
  libmkl_intel_lp64.a 
  libmkl_intel_ilp64.a 
  libmkl_intel_thread.a
  libmkl_core.a
)

set(MKL_FOUND "NO")
if(MKL_LIBRARY AND MKL_INCLUDE_DIR)
  set(MKL_FOUND "YES")
endif(MKL_LIBRARY AND MKL_INCLUDE_DIR)

if(MKL_FOUND)
  message(STATUS "Found intel mkl")
  message(STATUS " ${MKL_LIBRARY}")
  message(STATUS " ${MKL_INCLUDE_DIR}")
else(MKL_FOUND)
  message(STATUS "Could not find intel mkl")
endif(MKL_FOUND)
mark_as_advanced(MKL_LIBRARY  MKL_INCLUDE_DIR)
