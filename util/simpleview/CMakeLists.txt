#*******************************************************************************
# Simple View
# Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# - show results using gnuplot
#
# *******************************************************************************/

project(SIMPLE_VIEW)
cmake_minimum_required(VERSION 2.8)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/../gnuplot-iostream ")

# set general g++ flags
# note: we always need -O3 because minlin_host doesn't compile without it
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

add_executable(simpleview "simpleview.cpp")

# link external libraries	
target_link_libraries(simpleview "-lboost_iostreams -lboost_system -lboost_filesystem -lboost_program_options")

# set the name of output file
#set_source_files_properties(simpleview
#				COMPILE_FLAGS "")

set_target_properties(simpleview PROPERTIES
	OUTPUT_NAME simpleview
)

