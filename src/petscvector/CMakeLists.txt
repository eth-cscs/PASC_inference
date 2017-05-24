project(PASC_INFERENCE_LIBRARY_PETSC)
cmake_minimum_required(VERSION 2.8.11)

# compile cpp files with petsc
if(${USE_PETSC})
	file(GLOB_RECURSE SOURCES_PETSC "${PASCINFERENCE_SRC}/petscvector/*.cpp")
	set_source_files_properties(${SOURCES_PETSC} COMPILE_FLAGS "${FLAGS_DEF_D}")
	add_library(${PROJECT_NAME} OBJECT ${SOURCES_PETSC})
endif()


#project(PASC_INFERENCE_LIBRARY_PETSC_CUDA)

if(${USE_PETSC})
	# compile cuda cu files
	if(${USE_CUDA})
		file(GLOB_RECURSE SOURCES_PETSC_CUDA "${PASCINFERENCE_SRC}/petscvector/*.cu") # locate all cu files in source directory
		set_source_files_properties(${SOURCES_PETSC_CUDA} COMPILE_FLAGS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\"")
		set(PASC_INFERENCE_PETSC_CUDA_SOURCES ${SOURCES_PETSC_CUDA} PARENT_SCOPE)

#		cuda_add_library(CUDASSS ${SOURCES_PETSC_CUDA})
#		install(FILES "${CMAKE_BINARY_DIR}/lib/lib${PROJECT_NAME}_CUDA.a" DESTINATION lib) 	
#		cuda_compile("MYCUDA_O" ${SOURCES_PETSC_CUDA})
#		set(PASC_INFERENCE_LIBRARY_PETSC_CUDA_O ${PASCINFERENCE_LIBRARY_PETSC_CUDA_O} PARENT_SCOPE)

#		message(${PASC_INFERENCE_LIBRARY_PETSC_CUDA})
	endif()
endif()

