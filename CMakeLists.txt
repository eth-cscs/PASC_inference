#--------------------------------------------------
# MinLin CMAKE sample
#
# Requirements:
# - defined variable MKLROOT
# - minlin folder with source files in CMAKE_CURRENT_SOURCE_DIR
# - 
#
#--------------------------------------------------
project(QPOPT_MINLIN_SAMPLE)
cmake_minimum_required(VERSION 2.8)

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
set(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
include_directories(${MKL_INCLUDE_DIR})
link_directories($ENV{MKLROOT}/lib/intel64)

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# MINLIN: define variables for include directories
set(MINLIN_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
message(STATUS "MINLIN_INCLUDE = ${MINLIN_INCLUDE}")
set(MINLIN_HOST_DEFS
	THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
	THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
	__host__=\ 
	__device__=\ 
	USE_MINLIN)
set(MINLIN_DEVICE_DEFS # we use -D here because this isn't added automatically
	-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
	-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA
	-DUSE_GPU
	-DUSE_MINLIN)
include_directories(${MINLIN_INCLUDE})

# MINLIN HOST EXECUTABLE
add_executable(string_host string.cpp)
target_link_libraries(string_host mkl_core mkl_gnu_thread mkl_rt boost_program_options)
set_property(
	TARGET string_host
	APPEND
	PROPERTY INCLUDE_DIRECTORIES
	${MINLIN_INCLUDE})
set_target_properties(string_host PROPERTIES
	OUTPUT_NAME string_host
	COMPILE_DEFINITIONS "${MINLIN_HOST_DEFS}"
	)

add_executable(string_bound_host string_bound.cpp)
target_link_libraries(string_bound_host mkl_core mkl_gnu_thread mkl_rt boost_program_options)
set_property(
	TARGET string_bound_host
	APPEND
	PROPERTY INCLUDE_DIRECTORIES
	${MINLIN_INCLUDE})
set_target_properties(string_bound_host PROPERTIES
	OUTPUT_NAME string_bound_host
	COMPILE_DEFINITIONS "${MINLIN_HOST_DEFS}"
	)

add_executable(string_eq_host string_eq.cpp)
target_link_libraries(string_eq_host mkl_core mkl_gnu_thread mkl_rt boost_program_options)
set_property(
	TARGET string_eq_host
	APPEND
	PROPERTY INCLUDE_DIRECTORIES
	${MINLIN_INCLUDE})
set_target_properties(string_eq_host PROPERTIES
	OUTPUT_NAME string_eq_host
	COMPILE_DEFINITIONS "${MINLIN_HOST_DEFS}"
	)

# MINLIN DEVICE EXECUTABLE
include(FindCUDA)
set(CUDA_PROPAGATE_HOST_FLAGS off) # if flags are passed with -Xcompiler, they also affect NVCC which doesn't understand all g++ flags we use
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # without this, cc is used instead of CC and all include paths have to be specified manually
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
set(CUDA_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-vla ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}") # add flags specific to build type
string(REPLACE "-std=c++11" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++11 from options
cuda_include_directories(${MINLIN_INCLUDE})

# add all samples
cuda_add_executable(string_device string.cu
	OPTIONS "${MINLIN_DEVICE_DEFS} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
	DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
target_link_libraries(string_device mkl_core mkl_gnu_thread mkl_rt boost_program_options cublas)
set_target_properties(string_device PROPERTIES
	OUTPUT_NAME string_device
	)

cuda_add_executable(string_bound_device string_bound.cu
	OPTIONS "${MINLIN_DEVICE_DEFS} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
	DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
target_link_libraries(string_bound_device mkl_core mkl_gnu_thread mkl_rt boost_program_options cublas)
set_target_properties(string_bound_device PROPERTIES
	OUTPUT_NAME string_bound_device
	)

cuda_add_executable(string_eq_device string_eq.cu
	OPTIONS "${MINLIN_DEVICE_DEFS} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
	DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
target_link_libraries(string_eq_device mkl_core mkl_gnu_thread mkl_rt boost_program_options cublas)
set_target_properties(string_eq_device PROPERTIES
	OUTPUT_NAME string_eq_device
	)
