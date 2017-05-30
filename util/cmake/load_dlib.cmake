
if(${USE_DLIB})
	message(STATUS "${Blue}loading Dlib library${ColourReset}")

	set(DLIB_USE_MKL_FFT 0)
	set(DLIB_USE_CUDA 0)
	set(DLIB_LINK_WITH_SQLITE3 0)
	option(DLIB_USE_MKL_FFT OFF)
	option(DLIB_USE_USE_MKL_FFT OFF)

	# set the root to DLib library
	set(DLIB_ROOT "${PASCINFERENCE_ROOT}/util/dlib/")

	# the compilation of external library has to be solved separately
	# if the library is build as a part of different project
#	if(NOT ${PASCINFERENCE_IN_PROJECT_BUILD})
		# include some funny cmake functions from DLib
#	if(NOT EXISTS "${CMAKE_BINARY_DIR}/pascinference_build/dlib_build/libdlib.a")
#		include("${DLIB_ROOT}/dlib/cmake/")
#	endif()
#	endif()

	# append to flags definitions
	set(FLAGS_DEF "-USE_DLIB ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_DLIB ${FLAGS_DEF_D}")

	# try to find metis library
	find_library(DLIB_LIB_LOCATION "dlib")

	if(NOT DLIB_LIB_LOCATION)
		# add metis library for compilation
		if (NOT TARGET project_dlib)
			# if the library doesn't exist, then we will compile it from source code
			message(STATUS "\n${Blue}DLIB library not found, it will be compiled ${ColourReset}\n")

			# include cmake module to be able to use ExternalProject
			include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
		
			ExternalProject_Add(project_dlib
				SOURCE_DIR ${PASCINFERENCE_ROOT}/util/dlib/
#				GIT_SUBMODULES util/dlib
#				URL ${PASCINFERENCE_ROOT}/util/dlib/
				PREFIX ${CMAKE_BINARY_DIR}/dlib_build
				CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Release --config Release"
#				EXCLUDE_FROM_ALL TRUE
				STEP_TARGETS build
				INSTALL_COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/dlib_build/src/project_dlib-build/dlib/libdlib.a ${CMAKE_BINARY_DIR}/lib/libdlib.a
			)

			add_library(dlib SHARED IMPORTED)
		endif()
		
	else()
		# library found
		message(STATUS "${Blue}DLIB library found in: ${DLIB_LIB_LOCATION} ${ColourReset}")
	endif()

	# include header files
	link_directories("${PASCINFERENCE_ROOT}/util/dlib/")
	include_directories("${PASCINFERENCE_ROOT}/util/dlib/")

	# Now it is safe to include other dlib infrustucture - it won't build dlib again.


	# link library
	set(LIBRARIES_DEF "dlib;${LIBRARIES_DEF}")

#	set(LIBRARIES_DEF "dlib;pthread;X11;${LIBRARIES_DEF}")

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_DLIB)
	printinfo_onoff("USE_DLIB\t\t\t" "${USE_DLIB}")
	if(${USE_DLIB})
		printinfo(" - DLIB_ROOT\t\t\t" "${DLIB_ROOT}")
	endif()

endmacro()

