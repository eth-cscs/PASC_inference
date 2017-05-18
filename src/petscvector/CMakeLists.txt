project(PASC_INFERENCE_LIBRARY_PETSC)
cmake_minimum_required(VERSION 2.8.11)

# compile cpp files with petsc
if(${USE_PETSC})
	file(GLOB_RECURSE SOURCES_PETSC "${PASCINFERENCE_SRC}/petscvector/*.cpp")
	set_source_files_properties(${SOURCES_PETSC} COMPILE_FLAGS "${FLAGS_DEF_D}")
	add_library(${PROJECT_NAME} OBJECT ${SOURCES_PETSC})

	# compile cuda cu files
	if(${USE_CUDA})
		file(GLOB_RECURSE SOURCES_CUDA "${PASCINFERENCE_SRC}/petscvector/*.cu") # locate all cu files in source directory
		set_source_files_properties(${SOURCES_CUDA} COMPILE_FLAGS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\"")
		cuda_add_library(${PROJECT_NAME} OBJECT ${SOURCES_CUDA})
	endif()
endif()

