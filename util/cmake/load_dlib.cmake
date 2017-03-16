
if(${USE_DLIB})
	message(STATUS "${Blue}loading Dlib library${ColourReset}")

	# set the root to DLib library
	set(DLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../util/dlib/")

	# include some funny cmake functions from DLib
	include("${DLIB_ROOT}/dlib/cmake/")

	# append to flags definitions
	set(FLAGS_DEF "-USE_DLIB ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_DLIB ${FLAGS_DEF_D}")

	# link library
	set(LIBRARIES_DEF "dlib::dlib;${LIBRARIES_DEF}")
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_DLIB)
	printinfo_onoff("USE_DLIB\t\t\t" "${USE_DLIB}")
	if(${USE_DLIB})
		printinfo(" - DLIB_ROOT\t\t\t" "${DLIB_ROOT}")
	endif()

endmacro()

