# http://glaros.dtc.umn.edu/gkhome/metis/metis/download
# 
# compile using: 
# cd util/metis/build
# cmake -DSHARED=ON ..
# make
#

if(${USE_METIS})
	message(STATUS "${Yellow}loading METIS${ColourReset}")

	set(METIS_INCLUDE_DIR "${PASCINFERENCE_ROOT}/util/metis/include")

	if(${USE_CUDA})
		include_directories(${METIS_INCLUDE_DIR})
		cuda_include_directories(${METIS_INCLUDE_DIR})
	else()
		include_directories(${METIS_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_METIS ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_METIS ${FLAGS_DEF_D}")

	# try to find metis library
	find_library(METIS_LIB_LOCATION "metis")

	if(NOT METIS_LIB_LOCATION)
		# add metis library for compilation
		if (NOT TARGET project_metis)
			# if the library doesn't exist, then we will compile it from source code
			message(STATUS "\n${Blue}METIS library not found, it will be compiled ${ColourReset}\n")

			# include cmake module to be able to use ExternalProject
			include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
		
			ExternalProject_Add(project_metis
				SOURCE_DIR ${PASCINFERENCE_ROOT}/util/metis/
#				CMAKE_COMMAND cmake
#				GIT_SUBMODULES util/minlin
#				URL ${PASCINFERENCE_ROOT}/util/metis
				PREFIX ${CMAKE_BINARY_DIR}/metis_build
				CMAKE_ARGS "-DSHARED=ON shared=1"
				EXCLUDE_FROM_ALL TRUE
				STEP_TARGETS build
				INSTALL_COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/metis_build/src/project_metis-build/libmetis/libmetis.so ${CMAKE_BINARY_DIR}/lib/libmetis.so
			)

			add_library(metis SHARED IMPORTED)
		endif()
		
		# add dependency to main library
#		add_dependencies(pascinference project_metis)
		
#		add_subdirectory(${PASCINFERENCE_ROOT}/util/metis/ metis_build)
	else()
		# library found
		message(STATUS "${Blue}METIS library found in: ${METIS_LIB_LOCATION} ${ColourReset}")
	endif()

	# include compiled! METIS library
#	link_directories("${CMAKE_BINARY_DIR}/metis_build/src/project_metis-build/libmetis/")

#	find_library(METIS_LIB_LOCATION "metis" PATH "${PASCINFERENCE_ROOT}/util/metis/build/libmetis/")
#	find_library(METIS_LIB_LOCATION "metis" PATH "${PASCINFERENCE_ROOT}/util/metis/build/libmetis/")

	#set(METIS_LIB_LOCATION "${CMAKE_SOURCE_DIR}/../util/metis/build/libmetis/libmetis.so")
	
#	if(NOT METIS_LIB_LOCATION)
		# if the library doesn't exist, then give error
#		message(FATAL_ERROR "\n${Red}METIS library not found, did you forget to compile it? : ${METIS_LIB_LOCATION} ${ColourReset}\n")
#	else()
		# library found
#		message(STATUS "${Blue}METIS library found in: ${METIS_LIB_LOCATION} ${ColourReset}")

#	endif()

	# link library
	set(LIBRARIES_DEF "metis;${LIBRARIES_DEF}")

endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_METIS)
	printinfo_onoff("USE_METIS\t\t\t" "${USE_METIS}")
	if(${USE_METIS})
		printinfo(" - METIS_INCLUDE_DIR\t\t" "${METIS_INCLUDE_DIR}")
		printinfo(" - METIS_LIB_LOCATION\t\t" "${METIS_LIB_LOCATION}")
	endif()
endmacro()

