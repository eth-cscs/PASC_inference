# This is a CMake file meant to be included via include()
# It will trigger a compilation of pascinference *in the project* 
# including it

cmake_minimum_required(VERSION 2.8.12)

set(PASCINFERENCE_IN_PROJECT_BUILD true)

# Determine the path to pascinference.
if(NOT ${PASCINFERENCE_ROOT})
	set(PASCINFERENCE_ROOT ${CMAKE_CURRENT_LIST_DIR})
endif()

# Don't add dlib if it's already been added to the cmake project
if (NOT TARGET pascinference)
    add_subdirectory(${PASCINFERENCE_ROOT} pascinference_build)
endif()

