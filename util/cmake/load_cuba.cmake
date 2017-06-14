
if(${USE_CUBA})
	message(STATUS "${Yellow}loading CUBA${ColourReset}")

	set(CUBA_INCLUDE_DIR "${PASCINFERENCE_ROOT}/util/cuba")

	if(${USE_CUDA})
		include_directories(${CUBA_INCLUDE_DIR})
		cuda_include_directories(${CUBA_INCLUDE_DIR})
	else()
		include_directories(${CUBA_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_CUBA ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_CUBA ${FLAGS_DEF_D}")

	# try to find metis library
	find_library(CUBA_LIB_LOCATION "cuba")

	if(NOT CUBA_LIB_LOCATION)
		# add cuba library for compilation
		if (NOT TARGET project_cuba)
			# if the library doesn't exist, then we will compile it from source code
			message(STATUS "\n${Blue}CUBA library not found, it will be compiled ${ColourReset}\n")

			# include cmake module to be able to use ExternalProject
			include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
		
			ExternalProject_Add(project_cuba
				SOURCE_DIR ${PASCINFERENCE_ROOT}/util/cuba/
#				CMAKE_COMMAND cmake
#				GIT_SUBMODULES util/minlin
#				URL ${PASCINFERENCE_ROOT}/util/metis
				PREFIX ${CMAKE_BINARY_DIR}/cuba_build
#				CONFIGURE_COMMAND "./configure" 
				CONFIGURE_COMMAND "./configure" 
				BUILD_COMMAND "make"
				UPDATE_COMMAND ""
				BUILD_IN_SOURCE 0
				INSTALL_COMMAND ${CMAKE_COMMAND} -E create_symlink ${PASCINFERENCE_ROOT}/util/cuba/libcuba.a ${CMAKE_BINARY_DIR}/lib/libcuba.a

			)

			add_library(cuba SHARED IMPORTED)
		endif()
		
	else()
		# library found
		message(STATUS "${Blue}CUBA library found in: ${CUBA_LIB_LOCATION} ${ColourReset}")
	endif()

	# link library
	set(LIBRARIES_DEF "cuba;${LIBRARIES_DEF}")

endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_CUBA)
	printinfo_onoff("USE_CUBA\t\t\t" "${USE_CUBA}")
	if(${USE_CUBA})
		printinfo(" - CUBA_INCLUDE_DIR\t\t" "${CUBA_INCLUDE_DIR}")
		printinfo(" - CUBA_LIB_LOCATION\t\t" "${CUBA_LIB_LOCATION}")
	endif()
endmacro()

