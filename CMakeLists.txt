#*******************************************************************************
# PASC INFERENCE library
#
# compile shared library
#
#
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2017
# lukas.pospisil@usi.ch
#
# *******************************************************************************/
#
# cmake -DFIND_PETSC=ON -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx ..
#


project(PASC_INFERENCE_LIBRARY)
cmake_minimum_required(VERSION 2.8.11)

set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}")
set(PASCINFERENCE_INCLUDE "${PASCINFERENCE_ROOT}/include/")
set(PASCINFERENCE_SRC "${PASCINFERENCE_ROOT}/src/")
set(PASCINFERENCE_CMAKE "${PASCINFERENCE_ROOT}/util/cmake/")

# default settings
option(USE_BOOST "USE_BOOST" ON)
option(USE_CUDA "USE_CUDA" OFF)
option(USE_MINLIN "USE_MINLIN" OFF)
option(USE_MKL "USE_MKL" OFF)
option(USE_PETSC "USE_PETSC" ON)
option(FIND_PETSC "FIND_PETSC" ON)
option(USE_PERMON "USE_PERMON" OFF)
option(USE_METIS "USE_METIS" ON)
option(USE_CRAYPOWER "USE_CRAYPOWER" ON)
option(USE_DLIB "USE_DLIB" OFF)

# include cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_CMAKE}" ${CMAKE_MODULE_PATH})
include(pascinference) # include common cmake functions

# write something funny
message(STATUS "${Blue}compiling PASCInference library${ColourReset}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC") # we are compiling shared library

# include header files
include_directories("${PASCINFERENCE_INCLUDE}")

# add files to library
add_subdirectory("${PASCINFERENCE_SRC}/general/")
add_subdirectory("${PASCINFERENCE_SRC}/boost/")
add_subdirectory("${PASCINFERENCE_SRC}/cuda/")

# create list of targets
set(PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_GENERAL>)
if(${USE_BOOST})
	list(APPEND PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_BOOST>)
endif()
if(${USE_CUDA})
	list(APPEND PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_CUDA>)
endif()


# add library from created list
add_library(pascinference SHARED ${PASC_TARGETS})


# print information about settings
include(printsetting)


