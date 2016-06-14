if(${USE_EDFLIB})
	message(STATUS "${Yellow}loading EDFlib${ColourReset}")

	set(EDFLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../util/EDFlib)

	if(${USE_CUDA})
		include_directories(${EDFLIB_INCLUDE_DIR})
		cuda_include_directories(${EDFLIB_INCLUDE_DIR})
	else()
		include_directories(${EDFLIB_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_EDFLIB ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_EDFLIB ${FLAGS_DEF_D}")

	# check if library exists
	find_library(EDFLIB_LIB_LOCATION "edflib" ${MYLIBDIR})
	
	if(${EDFLIB_LIB_LOCATION})
		# library found
	
	else()
		# if the library doesn't exist, then compile it
		message(STATUS "${Blue}compiling EDFlib library${ColourReset}")

#		if(${USE_GPU})
#			message(STATUS "${Blue}- will be compiled using CUDA compiler${ColourReset}")
#			CUDA_ADD_LIBRARY(libpascinference ${PASCINFERENCE_SRC}/libpascinference.cu OPTIONS ${FLAGS_DEF_D} SHARED )
#		else()
			message(STATUS "${Blue}- will be compiled using C++ compiler${ColourReset}")
			add_library(edflib SHARED "${EDFLIB_INCLUDE_DIR}/edflib.c")
			set_source_files_properties("${EDFLIB_INCLUDE_DIR}/edflib.c"
				COMPILE_FLAGS "-Wall -Wextra -Wshadow -Wformat-nonliteral -Wformat-security -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE"
				PROPERTIES
				 LIBRARY_OUTPUT_DIRECTORY "${MYLIBDIR}")		
#		endif()

		# the library should be compiled first, then compile other stuff
		set(COMPILE_FIRST "edflib;${COMPILE_FIRST}")
	endif()

	# link library
	set(LIBRARIES_DEF "edflib;${LIBRARIES_DEF}")

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_EDFLIB)
	printinfo_onoff("USE_EDFLIB\t\t\t" "${USE_EDFLIB}")
	if(${USE_EDFLIB})
		printinfo(" - EDFLIB_INCLUDE_DIR\t\t" "${EDFLIB_INCLUDE_DIR}")
		printinfo(" - EDFLIB_LIB_LOCATION\t\t" "${EDFLIB_LIB_LOCATION}")
	endif()
endmacro()

