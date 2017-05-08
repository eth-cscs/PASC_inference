# This is a CMake file meant to be included via include()
# It will trigger a compilation of pascinference *in the project* 
# including it

cmake_minimum_required(VERSION 2.8.12)

set(PASCINFERENCE_IN_PROJECT_BUILD true)

#if (POLICY CMP0054)
#    cmake_policy(SET CMP0054 NEW)
#endif()

# Determine the path to pascinference.
if(NOT ${PASCINFERENCE_ROOT})
	set(PASCINFERENCE_ROOT ${CMAKE_CURRENT_LIST_DIR})
endif()

	#include(${dlib_path}/cmake_utils/add_global_compiler_switch.cmake)
#include(${dlib_path}/cmake_utils/use_cpp_11.cmake)


# This is really optional, but nice.  It will make sure the build mode
# created by cmake is always release by default.
#include(${dlib_path}/cmake_utils/release_build_by_default)


# Don't add dlib if it's already been added to the cmake project
if (NOT TARGET pascinference)
    add_subdirectory(${PASCINFERENCE_ROOT} pascinference_build)
endif()

