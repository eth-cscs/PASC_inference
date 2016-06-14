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
option(USE_EDFLIB "USE_EDFLIB" ON)
option(COMPILE_PASCINFERENCE "COMPILE_PASCINFERENCE" OFF)

# include cmake functions
include(common) # include common

# print information about settings
include(printsetting)

# prepare directories
make_directory("results")

# decide which example to compile
option(TSSOLVER_GLOBAL_KMEANS "TSSOLVER_GLOBAL_KMEANS" OFF)
option(TSSOLVER_GLOBAL_VARX "TSSOLVER_GLOBAL_VARX" OFF)
option(PROJECTION_TEST "PROJECTION_TEST" OFF)
option(EDF "EDF" OFF)

message(STATUS "${Yellow}TSSOLVER_GLOBAL_KMEANS : ${Green}${TSSOLVER_GLOBAL_KMEANS}${Yellow}${ColourReset}")
message(STATUS "${Yellow}TSSOLVER_GLOBAL_VARX   : ${Green}${TSSOLVER_GLOBAL_VARX}${Yellow}${ColourReset}")
message(STATUS "${Yellow}PROJECTION_TEST        : ${Green}${PROJECTION_TEST}${Yellow}${ColourReset}")
message(STATUS "${Yellow}EDF                    : ${Green}${EDF}${Yellow}${ColourReset}")

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

if(${TSSOLVER_GLOBAL_VARX})
	# this is VARX global test
	if(${USE_CUDA})
		pascadd_executable("tssolver_global_varx.cu" "tssolver_global_varx")
	else()
		pascadd_executable("tssolver_global_varx.cpp" "tssolver_global_varx")
	endif()
endif()

if(${PROJECTION_TEST})
	# this is projection test
	if(${USE_CUDA})
		pascadd_executable("projection_test.cu" "projection_test")
	else()
		pascadd_executable("projection_test.cpp" "projection_test")
	endif()
endif()

if(${EDF})
	# this is edf test
	if(${USE_CUDA})
		pascadd_executable("edf.cu" "edf")
	else()
		pascadd_executable("edf.cpp" "edf")
	endif()
endif()


