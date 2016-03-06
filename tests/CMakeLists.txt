#*******************************************************************************
# PASC INFERENCE library
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# *******************************************************************************/

project(PASC_INFERENCE)
cmake_minimum_required(VERSION 2.8)

# add cmake functions
set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}/../")

# add cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/" ${CMAKE_MODULE_PATH})

include(common) # include common
include(load_gpu) # GPU
include(load_petsc) # PETSC
include(load_petscvector) # PetscVector
include(load_mkl) # MKL
include(load_minlin) # MinLin
include(load_pascinference) # PascInference


# print information about settings
include(printsetting)


















