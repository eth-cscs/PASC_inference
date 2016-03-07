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
option(TEST_GIVENMATRIX "TEST_GIVENMATRIX" OFF)
option(TEST_LAPLACEMATRIX "TEST_LAPLACEMATRIX" OFF)
option(TEST_CG "TEST_CG" OFF)
option(TEST_SPGQP "TEST_SPGQP" OFF)
option(TEST_KMEANS "TEST_KMEANS" ON)


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
pascadd_executable("test_algebra.cu" "test_algebra")













