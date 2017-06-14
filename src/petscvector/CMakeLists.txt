project(PASC_INFERENCE_LIBRARY_PETSC)
cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/")

# compile cpp files with petsc
if(${USE_PETSC})
	file(GLOB_RECURSE SOURCES_PETSC "${PASCINFERENCE_SRC}/petscvector/*.cpp")
#	file(GLOB_RECURSE SOURCES_PETSC "${PASCINFERENCE_SRC}/petscvector/*initialize.cpp")

	set_source_files_properties(${SOURCES_PETSC} COMPILE_FLAGS "${FLAGS_DEF_D}")
	add_library(pascinference_petsc SHARED ${SOURCES_PETSC})

	# link external libraries	
#	target_link_libraries(pascinference_petsc ${PETSC_LIBRARIES})

	# compile cuda cu files
	if(${USE_CUDA})
#		file(GLOB_RECURSE SOURCES_PETSC_CUDA "${PASCINFERENCE_SRC}/petscvector/*initialize.cu") # locate all cu files in source directory
		file(GLOB_RECURSE SOURCES_PETSC_CUDA "${PASCINFERENCE_SRC}/petscvector/*.cu")

#		set_source_files_properties(${SOURCES_PETSC_CUDA} COMPILE_FLAGS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\"")
		cuda_add_library(pascinference_petsc_cuda SHARED ${SOURCES_PETSC_CUDA} 
			OPTIONS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"-fPIC ${CUDA_CXX_FLAGS} \""
			DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

		# add dependency to original petsc library
#		add_dependencies(pascinference_petsc pascinference_petsc_cuda)

		# link external libraries	
#		target_link_libraries(pascinference_petsc_cuda ${LIBRARIES_DEF})
	endif()

	# if there is cuba, then add it as a dependency
	if(${USE_CUBA} AND TARGET project_cuba)
		add_dependencies(pascinference_petsc project_cuba)
	endif()

	# if there is dlib, then add it as a dependency
	if(${USE_DLIB} AND TARGET project_dlib)
		add_dependencies(pascinference_petsc project_dlib)
	endif()
		
	# if there is metis, then add it as a dependency
	if(${USE_METIS} AND TARGET project_metis)
		add_dependencies(pascinference_petsc project_metis)
	endif()

	# if there is permon, then add it as a dependency
	if(${USE_PERMON} AND TARGET project_permon)
		add_dependencies(pascinference_petsc project_permon)
	endif()

endif()

