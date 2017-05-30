# PERMON is a set of solvers combining quadratic programming and domain decomposition methods. It makes use of and extends the PETSc framework for numerical computations.
# http://permon.it4i.cz/
# https://github.com/haplav/permon
# 
# compile using:
# export PERMON_DIR=[PASC_INFERENCE_DIR]/util/permon 
# cd $PERMON_DIR
# make
#


if(${USE_PERMON})
	message(STATUS "${Yellow}loading PERMON${ColourReset}")

	set(FLAGS_DEF "-USE_PERMON ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_PERMON ${FLAGS_DEF_D}")

	# try to find PERMON library
	find_library(PERMON_LIB_LOCATION "permon" PATH "$ENV{PERMON_DIR}/$ENV{PETSC_ARCH}/lib/")

	if(NOT PERMON_LIB_LOCATION)
		# add metis library for compilation
		if (NOT TARGET project_permon)
			# if the library doesn't exist, then we will compile it from source code
			message(STATUS "\n${Blue}PERMON library not found, it will be compiled ${ColourReset}\n")

			# include cmake module to be able to use ExternalProject
			include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
		
			set( ENV{PERMON_DIR} "${PASCINFERENCE_ROOT}/util/permon")
		
			ExternalProject_Add(project_permon
				SOURCE_DIR ${PASCINFERENCE_ROOT}/util/permon/
#				GIT_SUBMODULES util/permon
#				URL ${PASCINFERENCE_ROOT}/util/permon
				PREFIX ${CMAKE_BINARY_DIR}/permon_build
#				CMAKE_ARGS "-DSHARED=ON shared=1 -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}/lib/"
#				EXCLUDE_FROM_ALL TRUE
#				STEP_TARGETS build
				CONFIGURE_COMMAND "" 
				BUILD_COMMAND make PERMON_DIR=$ENV{PERMON_DIR}
				UPDATE_COMMAND ""
				BUILD_IN_SOURCE 1
				INSTALL_COMMAND ${CMAKE_COMMAND} -E create_symlink ${PASCINFERENCE_ROOT}/util/permon/$ENV{PETSC_ARCH}/lib/libpermon.so ${CMAKE_BINARY_DIR}/lib/libpermon.so
			)
		
			add_library(permon SHARED IMPORTED)
			
		endif()

	else()
		# library found
		message(STATUS "${Blue}PERMON library found in: ${PERMON_LIB_LOCATION} ${ColourReset}")
	endif()

#	link_directories("$ENV{PERMON_DIR}/$ENV{PETSC_ARCH}/lib/")

	# include header files
	set(PERMON_INCLUDE_DIR "$ENV{PERMON_DIR}/include")

	if(${USE_CUDA})
		include_directories(${PERMON_INCLUDE_DIR})
		cuda_include_directories(${PERMON_INCLUDE_DIR})
	else()
		include_directories(${PERMON_INCLUDE_DIR})
	endif()

	# link library
	set(LIBRARIES_DEF "permon;${LIBRARIES_DEF}")

endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PERMON)
	printinfo_onoff("USE_PERMON\t\t\t" "${USE_PERMON}")
	if(${USE_PERMON})
		printinfo(" - PERMON_DIR\t\t\t" "$ENV{PERMON_DIR}")
		printinfo(" - PERMON_INCLUDE_DIR\t\t" "${PERMON_INCLUDE_DIR}")
		printinfo(" - PERMON_LIB_LOCATION\t\t" "${PERMON_LIB_LOCATION}")
	endif()
endmacro()

