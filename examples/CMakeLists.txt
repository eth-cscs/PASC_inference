#*******************************************************************************
# PASC INFERENCE library examples
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# *******************************************************************************/
#
# cmake -DFIND_PETSC=ON -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx ..
#

project(PASC_INFERENCE)
cmake_minimum_required(VERSION 2.8)

#decide which example to compile
option(TSSOLVER_GLOBAL_KMEANS "TSSOLVER_GLOBAL_KMEANS" OFF)
option(TSSOLVER_GLOBAL_KMEANS_LARGE "TSSOLVER_GLOBAL_KMEANS_LARGE" ON)
option(TSSOLVER_GLOBAL_VARX "TSSOLVER_GLOBAL_VARX" OFF)

# include additional example headers
include_directories("include")

# set the root to the library
set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}/../")

# include cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/" ${CMAKE_MODULE_PATH})

#default settings
option(USE_BOOST "USE_BOOST" ON)
option(USE_CUDA "USE_CUDA" OFF)
option(USE_GPU "USE_GPU" OFF)
option(USE_MINLIN "USE_MINLIN" OFF)
option(USE_MKL "USE_MKL" OFF)
option(USE_PETSC "USE_PETSC" ON)
option(FIND_PETSC "FIND_PETSC" OFF)
option(USE_PETSCVECTOR "USE_PETSCVECTOR" ON)
option(COMPILE_PASCINFERENCE "COMPILE_PASCINFERENCE" OFF)

# include cmake functions
include(common) # include common
include(load_cuda) # CUDA
include(load_gpu) # GPU
include(load_petsc) # PETSC
include(load_petscvector) # PetscVector
include(load_boost) # BOOST
include(load_mkl) # MKL
include(load_minlin) # MinLin
include(load_pascinference) # PascInference

# print information about settings
include(printsetting)

# prepare directories
make_directory("results")

# add example executable files
if(${TSSOLVER_GLOBAL_KMEANS})
	# this is VARX global test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("tssolver_global_kmeans.cu" "tssolver_global_kmeans")
	else()
		# otherwise compile as .cpp
		pascadd_executable("tssolver_global_kmeans.cpp" "tssolver_global_kmeans")
	endif()
endif()

if(${TSSOLVER_GLOBAL_KMEANS_LARGE})
	# this is VARX global test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("tssolver_global_kmeans_large.cu" "tssolver_global_kmeans_large")
	else()
		# otherwise compile as .cpp
		pascadd_executable("tssolver_global_kmeans_large.cpp" "tssolver_global_kmeans_large")
	endif()
endif()

if(${TSSOLVER_GLOBAL_VARX})
	# this is VARX global test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("tssolver_global_varx.cu" "tssolver_global_varx")
	else()
		# otherwise compile as .cpp
		pascadd_executable("tssolver_global_varx.cpp" "tssolver_global_varx")
	endif()
endif()

