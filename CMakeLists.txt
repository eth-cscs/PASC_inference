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
include(common) # include common cmake functions
include(load_boost) # BOOST
include(load_cuda) # CUDA
include(load_petsc) # PETSC
include(load_permon) # PERMON
include(load_mkl) # MKL
include(load_minlin) # MinLin
include(load_metis) # METIS
include(load_craypower) # for measuring power consumption on Piz Daint
include(load_dlib) # DLIB

message(STATUS "${Blue}compiling PASCInference library${ColourReset}")

# include header files
include_directories("${PASCINFERENCE_INCLUDE}")

# compile general functions
file(GLOB LIBPASCINFERENCE_SOURCES "${PASCINFERENCE_SRC}/general/*.cpp") # locate all cpp files in source directory
set_source_files_properties(${LIBPASCINFERENCE_SOURCES}	COMPILE_FLAGS ${FLAGS_DEF_D})
add_library(pascinference_general SHARED ${LIBPASCINFERENCE_SOURCES})

# compile cuda functions
if(${USE_CUDA})
	file(GLOB LIBPASCINFERENCE_CUDA_SOURCES "${PASCINFERENCE_SRC}/cuda/*.cu")
	set_source_files_properties(${LIBPASCINFERENCE_CUDA_SOURCES} COMPILE_FLAGS ${FLAGS_DEF_D})
	cuda_add_library(pascinference_cuda SHARED ${LIBPASCINFERENCE_CUDA_SOURCES})
endif()

# compile petsc functions
if(${USE_PETSC})
#	file(GLOB LIBPASCINFERENCE_PETSC_SOURCES "${PASCINFERENCE_SRC}/petsc/*.cpp")
#	set_source_files_properties(${LIBPASCINFERENCE_PETSC_SOURCES} COMPILE_FLAGS ${FLAGS_DEF_D})
#	add_library(pascinference_petsc SHARED ${LIBPASCINFERENCE_PETSC_SOURCES})
endif()



# print information about settings
include(printsetting)


