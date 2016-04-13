#*******************************************************************************
# PASC INFERENCE library examples
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# *******************************************************************************/

project(PASC_INFERENCE)
cmake_minimum_required(VERSION 2.8)

#decide which example to compile
option(TEST_ALGEBRA "TEST_ALGEBRA" OFF)
option(TEST_CG "TEST_CG" OFF)
option(TEST_FILECRSMATRIX "TEST_FILECRSMATRIX" OFF)
option(TEST_LAPLACEMATRIX "TEST_LAPLACEMATRIX" OFF)
option(TEST_SPGQP "TEST_SPGQP" OFF)
option(TEST_KMEANS "TEST_KMEANS" OFF)
option(TEST_VARX "TEST_VARX" OFF)
option(BENCHMARK1 "BENCHMARK1" OFF)
option(BENCHMARK2 "BENCHMARK2" ON)

# include additional example headers
include_directories("include")

# set the root to the library
set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}/../")

# include cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/" ${CMAKE_MODULE_PATH})

include(common) # include common
include(load_cuda) # CUDA
include(load_gpu) # GPU
include(load_petsc) # PETSC
include(load_petscvector) # PetscVector
include(load_mkl) # MKL
include(load_minlin) # MinLin
include(load_pascinference) # PascInference

# print information about settings
include(printsetting)


# add example executable files
if(${TEST_ALGEBRA})
	pascadd_executable("test_algebra.cu" "test_algebra")
endif()

if(${TEST_FILECRSMATRIX})
	pascadd_executable("test_filecrsmatrix.cu" "test_filecrsmatrix")
endif()

if(${TEST_LAPLACEMATRIX})
	pascadd_executable("test_laplacematrix.cu" "test_laplacematrix")
endif()

if(${TEST_CG})
	pascadd_executable("test_cg.cu" "test_cg")
endif()

if(${TEST_SPGQP})
	pascadd_executable("test_spgqp.cu" "test_spgqp")
endif()

if(${TEST_KMEANS})
	# this is KMEANS test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_kmeans.cu" "test_kmeans")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_kmeans.cpp" "test_kmeans")
	endif()
endif()

if(${TEST_VARX})
	# this is VARX test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_varx.cu" "test_varx")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_varx.cpp" "test_varx")
	endif()
endif()

if(${BENCHMARK1})
	# this is VARX test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("benchmark1.cu" "benchmark1")
	else()
		# otherwise compile as .cpp
		pascadd_executable("benchmark1.cpp" "benchmark1")
	endif()
endif()

if(${BENCHMARK2})
	# this is VARX test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("benchmark2.cu" "benchmark2")
	else()
		# otherwise compile as .cpp
		pascadd_executable("benchmark2.cpp" "benchmark2")
	endif()
endif()



