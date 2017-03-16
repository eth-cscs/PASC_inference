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

# include additional example headers
include_directories("include")

# set the root to the library
set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}/../")

# include cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/" ${CMAKE_MODULE_PATH})

# default settings
option(USE_BOOST "USE_BOOST" ON)
option(USE_CUDA "USE_CUDA" OFF)
option(USE_GPU "USE_GPU" OFF)
option(USE_MINLIN "USE_MINLIN" OFF)
option(USE_MKL "USE_MKL" OFF)
option(USE_PETSC "USE_PETSC" ON)
option(FIND_PETSC "FIND_PETSC" OFF)
option(USE_PETSCVECTOR "USE_PETSCVECTOR" ON)
option(USE_PERMON "USE_PERMON" OFF)
option(USE_METIS "USE_METIS" ON)
option(USE_CRAYPOWER "USE_CRAYPOWER" ON)
option(USE_DLIB "USE_DLIB" OFF)
option(COMPILE_PASCINFERENCE "COMPILE_PASCINFERENCE" OFF)

# include cmake functions
include(common) # include common

# print information about settings
include(printsetting)

# prepare directories
make_directory("results")
make_directory("data")
make_directory("scripts")
make_directory("batch")
make_directory("batch_out")
make_directory("shortinfo")
make_directory("log")

# --------- include test cmake functions --------
include("util/util.cmake")
include("test_signal1D/test_signal1D.cmake")
include("test_entropy/test_entropy.cmake")
include("test_image/test_image.cmake")
include("test_neuro/test_neuro.cmake")
include("test_classes/test_classes.cmake")
message("\n----------------------------------------------------------\n")

