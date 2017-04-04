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

project(PASC_INFERENCE_TESTS)
cmake_minimum_required(VERSION 2.8.11)

# set PASC_Inference library path
set(PASCINFERENCE_ROOT "${CMAKE_SOURCE_DIR}/../")

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

# include pascinference cmake functions
set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/" ${CMAKE_MODULE_PATH})
include(pascinference) # include common cmake functions

# include PASC_Inference header files
include_directories("${PASCINFERENCE_ROOT}/include")

# prepare directories
make_directory("results")
make_directory("data")
make_directory("scripts")
make_directory("batch")
make_directory("batch_out")
make_directory("shortinfo")
make_directory("log")

# prepare macros for simple example compilation
macro(TESTADD_EXECUTABLE filename outname)
	# choose compiler subject to file extension
	get_filename_component(FILE_EXT ${filename} EXT)
	
	if(${FILE_EXT} MATCHES ".cu")
		# --- compile with CUDA ---
		if(NOT ${USE_CUDA})
			message(FATAL_ERROR "${Red}Cannot compile .cu file without USE_CUDA=ON!${ColourReset}")
		endif()
		
		# add executable file
		cuda_add_executable(${outname} ${filename}
			OPTIONS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\""
			DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	
		# set the name of output file
		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
		)
		
	endif()	

	if(${FILE_EXT} MATCHES ".cpp")
		# --- compile with g++ ---

		# add executable file
		add_executable(${outname} ${filename})

		set_source_files_properties(${filename}
				COMPILE_FLAGS "${FLAGS_DEF_D}")

		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
#			DEBUG ${CMAKE_CXX_FLAGS_DEBUG}
		)
	endif()	

	# link pasc inference library
	target_link_libraries(${outname} ${LIBRARIES_DEF})

endmacro()


# --------- include test cmake functions --------
include("test_util/test_util.cmake")
include("test_signal1D/test_signal1D.cmake")
#include("test_entropy/test_entropy.cmake")
#include("test_image/test_image.cmake")
#include("test_neuro/test_neuro.cmake")
#include("test_classes/test_classes.cmake")
message("\n----------------------------------------------------------\n")

