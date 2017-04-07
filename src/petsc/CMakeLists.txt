project(PASC_INFERENCE_LIBRARY_PETSC)
cmake_minimum_required(VERSION 2.8.11)

# compile cpp files with petsc
if(${USE_PETSC})
	file(GLOB_RECURSE SOURCES_PETSC "${PASCINFERENCE_SRC}/petsc/*.cpp")
	set_source_files_properties(${SOURCES_PETSC} COMPILE_FLAGS "${FLAGS_DEF_D}")
	add_library(${PROJECT_NAME} OBJECT ${SOURCES_PETSC})
endif()

