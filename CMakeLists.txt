cmake_minimum_required(VERSION 2.8)
project(MEX)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${PROJECT_SOURCE_DIR}/cmake-compat)

if(CMAKE_VERSION VERSION_LESS "3.3") 
  message(WARNING "CMake version less than 3.3\n"
    "using local FindMatlab module extracted from v3.3")
  # Where we keep FindMatlab module
  set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake-compat)
endif()

find_package(Matlab REQUIRED MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  message(FATAL_ERROR "MATLAB not found.")
else()
  message("MATLAB found.")
endif()
include_directories(
  ${Matlab_INCLUDE_DIRS}
)

find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
  message(FATAL_ERROR "GSL not found.")
else()
  message("GSL found.")
endif()

# Assuming gcc for now
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  # -ffast-math is for some reason not working with gcc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -O3 -march=native -fopenmp")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  link_directories(/NOBACKUP/davoudss/intel/lib/intel64/)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -fast -O3 -march=native -openmp")
endif()

find_package(MKL)
if(MKL_FOUND)
  include_directories(
    ${MKL_INCLUDE_DIR}
    )

  link_directories(
    ${MKL_LIBRARY}
    )
endif()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # Assuming gcc for now
  # -ffast-math is for some reason not working with gcc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -O3 -march=native -fopenmp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -O3 -ipo -fast -std=c99 -march=native -openmp")
endif()

# Add modules with MEX to be built
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_direct")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_Stokes_direct/mex")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_fast_gridding")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_Stokes/mex")  
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_Stresslet/mex")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_Rotlet/mex")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE2P_Stokes/mex")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE1P/mex")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE0P/mex")

# Enable running Matlab test suite through ctest
enable_testing()
matlab_add_unit_test(
  NAME matlab_test_suite
  UNITTEST_FILE "${PROJECT_SOURCE_DIR}/run_unit_tests.m"
  UNITTEST_PRECOMMAND "cd ${PROJECT_SOURCE_DIR}"
  NO_UNITTEST_FRAMEWORK
)
