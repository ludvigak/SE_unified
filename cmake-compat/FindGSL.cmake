# Try to find gnu scientific library GSL
# See 
# http://www.gnu.org/software/gsl/  and 
# http://gnuwin32.sourceforge.net/packages/gsl.htm
#
# Based on a script of Felix Woelk and Jan Woetzel
# (www.mip.informatik.uni-kiel.de)
# 
# It defines the following variables:
#  GSL_FOUND - system has GSL lib
#  GSL_INCLUDE_DIRS - where to find headers 
#  GSL_LIBRARIES - full path to the libraries
#  GSL_LIBRARY_DIRS, the directory where the PLplot library is found.

#  CMAKE_GSL_CXX_FLAGS  = Unix compiler flags for GSL, essentially "`gsl-config --cxxflags`"
#  GSL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
#  GSL_EXE_LINKER_FLAGS = rpath on Unix

SET(GSL_CONFIG_PREFER_PATH
  "$ENV{GSL_DIR}/bin"
  "$ENV{GSL_DIR}"
  "$ENV{GSL_HOME}/bin"
  "$ENV{GSL_HOME}"
  CACHE STRING "preferred path to GSL (gsl-config)")
FIND_PROGRAM(GSL_CONFIG gsl-config
  ${GSL_CONFIG_PREFER_PATH}
  /usr/bin/
  )
# MESSAGE("DBG GSL_CONFIG ${GSL_CONFIG}")

IF (GSL_CONFIG)
  # set CXXFLAGS to be fed into CXX_FLAGS by the user:
  SET(GSL_CXX_FLAGS "`${GSL_CONFIG} --cflags`")

  # set INCLUDE_DIRS to prefix+include
  EXEC_PROGRAM(${GSL_CONFIG}
    ARGS --prefix
    OUTPUT_VARIABLE GSL_PREFIX)
  SET(GSL_INCLUDE_DIR ${GSL_PREFIX}/include CACHE STRING INTERNAL)

  # set link libraries and link flags
  EXEC_PROGRAM(${GSL_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE GSL_LIBRARIES)
  SET(GSL_LIBRARIES "`${GSL_CONFIG} --libs`")

  # extract link dirs for rpath
  EXEC_PROGRAM(${GSL_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE GSL_CONFIG_LIBS )

  # split off the link dirs (for rpath)
  # use regular expression to match wildcard equivalent "-L*<endchar>"
  # with <endchar> is a space or a semicolon
  STRING(REGEX MATCHALL "[-][L]([^ ;])+"
    GSL_LINK_DIRECTORIES_WITH_PREFIX
    "${GSL_CONFIG_LIBS}" )
  #      MESSAGE("DBG  GSL_LINK_DIRECTORIES_WITH_PREFIX=${GSL_LINK_DIRECTORIES_WITH_PREFIX}")

  # remove prefix -L because we need the pure directory for LINK_DIRECTORIES

  IF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
    STRING(REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES ${GSL_LINK_DIRECTORIES_WITH_PREFIX} )
  ENDIF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
  SET(GSL_EXE_LINKER_FLAGS "-Wl,-rpath,${GSL_LINK_DIRECTORIES}" CACHE STRING INTERNAL)
  #      MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LINK_DIRECTORIES}")
  #      MESSAGE("DBG  GSL_EXE_LINKER_FLAGS=${GSL_EXE_LINKER_FLAGS}")

  #      ADD_DEFINITIONS("-DHAVE_GSL")
  #      SET(GSL_DEFINITIONS "-DHAVE_GSL")
  MARK_AS_ADVANCED(
    GSL_CXX_FLAGS
    GSL_INCLUDE_DIR
    GSL_LIBRARIES
    GSL_LINK_DIRECTORIES
    GSL_DEFINITIONS
    )
  #MESSAGE(STATUS "Using GSL from ${GSL_PREFIX}")

ELSE(GSL_CONFIG)
  MESSAGE("FindGSL.cmake: gsl-config not found. Please set it manually. GSL_CONFIG=${GSL_CONFIG}")
ENDIF(GSL_CONFIG)

IF(GSL_LIBRARIES)
  IF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
    SET(GSL_FOUND 1)
  ENDIF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
ELSE(GSL_LIBRARIES)
  IF (GSL_FIND_REQUIRED)
    message(SEND_ERROR "FindGSL.cmake: Unable to find the required GSL libraries")
  ENDIF(GSL_FIND_REQUIRED)
ENDIF(GSL_LIBRARIES)