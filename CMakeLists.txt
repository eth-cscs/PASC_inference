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

if(NOT ${PASCINFERENCE_ROOT} AND NOT ${PASCINFERENCE_IN_PROJECT_BUILD})
	set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}")
endif()
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

# print info about settings
if(NOT ${PASCINFERENCE_IN_PROJECT_BUILD})
	include(printsetting)
endif()

# write something funny
message(STATUS "${Blue}compiling PASCInference library${ColourReset}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC") # we are compiling shared library

# include header files
include_directories("${PASCINFERENCE_INCLUDE}")

# create list of general targets
set(PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_GENERAL>)
add_subdirectory("${PASCINFERENCE_SRC}/general/" ${CMAKE_BINARY_DIR}/pascinference/general)

# add seqarrayvector
add_subdirectory("${PASCINFERENCE_SRC}/seqarrayvector/" ${CMAKE_BINARY_DIR}/pascinference/seqarrayvector)
list(APPEND PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_SEQARRAYVECTOR>)

# add other external stuff
if(${USE_BOOST})
	add_subdirectory("${PASCINFERENCE_SRC}/boost/" ${CMAKE_BINARY_DIR}/pascinference/boost)
	list(APPEND PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_BOOST>)
endif()
if(${USE_PETSC})
	add_subdirectory("${PASCINFERENCE_SRC}/petscvector/" ${CMAKE_BINARY_DIR}/pascinference/petscvector)
	list(APPEND PASC_TARGETS $<TARGET_OBJECTS:PASC_INFERENCE_LIBRARY_PETSC>)
	if(${USE_CUDA})
		set(PASC_TARGETS_CUDA ${PASC_INFERENCE_PETSC_CUDA_SOURCES})
#		message("kunda: ${SOURCES_PETSC_CUDA}")
#		list(APPEND PASC_TARGETS ${MYCUDA_O})
#		target_link_libraries(pascinference CUDASSS)
	endif()
endif()


# add library from created list
add_library(pascinference SHARED ${PASC_TARGETS})
if(${USE_CUDA})
	cuda_add_library(pascinference_cuda SHARED ${PASC_TARGETS_CUDA})
endif()
#install(TARGETS pascinference DESTINATION lib)

