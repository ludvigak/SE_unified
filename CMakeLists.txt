cmake_minimum_required(VERSION 2.8)
project(MEX)

if(CMAKE_VERSION VERSION_LESS "3.3") 
  message(WARNING "CMake version less than 3.3\n"
    "using local FindMatlab module extracted from v3.3")
  # Where we keep FindMatlab module
  set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake)
endif()

find_package(Matlab REQUIRED)
if(NOT Matlab_FOUND)
  message(FATAL_ERROR "MATLAB not found.")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

# Assuming gcc for now
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -O3 -ffast-math -march=native -fopenmp")

include_directories(
  ${Matlab_INCLUDE_DIRS}
)

add_subdirectory("${PROJECT_SOURCE_DIR}/SE_fast_gridding")
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_Stresslet/mex")
