#*******************************************************************************
# PASC INFERENCE library
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# - uses MINLIN library with QPOPT solvers
#
# *******************************************************************************/

project(PASC_INFERENCE)
cmake_minimum_required(VERSION 2.8)

#load options from command line
option(USE_GPU "USE_GPU" ON)

#define some colors for fun
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
endif()

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

# include directories to this project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# set general g++ flags
# note: we always need -O3 because minlin_host doesn't compile without it
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")

# define variables for mkl include directories
if(NOT DEFINED ENV{MKLROOT})
  message(FATAL_ERROR "${Red}MKLROOT has to be specified!${ColourReset}")
  return()
endif()
set(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
message(STATUS "${Yellow}MKL_INCLUDE = ${MKL_INCLUDE_DIR}${ColourReset}")
include_directories(${MKL_INCLUDE_DIR})
link_directories($ENV{MKLROOT}/lib/intel64)

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# MINLIN: define variables for include directories
set(MINLIN_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/minlin/include)
message(STATUS "${Yellow}MINLIN_INCLUDE = ${MINLIN_INCLUDE}${ColourReset}")
set(MINLIN_HOST_DEFS
	THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
	THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
	__host__=\ 
	__device__=\ 
	USE_MINLIN)
set(MINLIN_DEVICE_DEFS # we use -D here because this isn't added automatically
	-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
	-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA
	-DUSE_MINLIN)
include_directories(${MINLIN_INCLUDE})

# MINLIN HOST EXECUTABLE
#add_executable(pascinference_host pascinference.cpp)
#target_link_libraries(pascinference_host mkl_core mkl_gnu_thread mkl_rt boost_program_options)
#set_property(
#	TARGET pascinference_host
#	APPEND
#	PROPERTY INCLUDE_DIRECTORIES
#	${MINLIN_INCLUDE})
#set_target_properties(pascinference_host PROPERTIES
#	OUTPUT_NAME pascinference_host
#	COMPILE_DEFINITIONS "${MINLIN_HOST_DEFS}"
#	)

# MINLIN DEVICE EXECUTABLE
include(FindCUDA)
set(CUDA_PROPAGATE_HOST_FLAGS off) # if flags are passed with -Xcompiler, they also affect NVCC which doesn't understand all g++ flags we use
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # without this, cc is used instead of CC and all include paths have to be specified manually
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
set(CUDA_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-vla ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}") # add flags specific to build type
string(REPLACE "-std=c++11" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++11 from options
cuda_include_directories(${MINLIN_INCLUDE})

# from source files create shared library
CUDA_ADD_LIBRARY(externals_common src/common.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_problem src/problem.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_data src/data.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_gamma src/gamma.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_theta src/theta.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_savevtk src/savevtk.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_qpsolver src/qpsolver.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_operations src/operations.cu
				 SHARED)
CUDA_ADD_LIBRARY(externals_projection src/projection.cu
				 SHARED)


# add main executable file
cuda_add_executable(pascinference_device pascinference.cu
	OPTIONS "${MINLIN_DEVICE_DEFS} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
	DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
target_link_libraries(pascinference_device mkl_core mkl_gnu_thread mkl_rt boost_program_options cublas)
set_target_properties(pascinference_device PROPERTIES
	OUTPUT_NAME pascinference
	)
target_link_libraries(pascinference_device 
	externals_common
	externals_problem
	externals_data
	externals_gamma
	externals_theta
	externals_savevtk
	externals_qpsolver
	externals_operations
	externals_projection
)

# add USE_GPU flag to all libs and exec files
if(${USE_GPU})
	message(STATUS "${Red}The code is compiled to run on GPU${ColourReset}")
	set(MINLIN_DEVICE_DEFS
		-DUSE_GPU)

else(${USE_GPU})
	message(STATUS "${Red}The code is compiled to run on CPU${ColourReset}")
endif(${USE_GPU})
	
# make output directory, here the results will be stored
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)	
