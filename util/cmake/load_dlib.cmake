
if(${USE_DLIB})
	message(STATUS "${Blue}loading Dlib library${ColourReset}")

	# disable MKL stuff if we are not using it
#	if(NOT ${USE_MKL})
		set(DLIB_USE_MKL_FFT 0)
#	endif()

#	if(NOT ${USE_CUDA})
		set(DLIB_USE_CUDA 0)
#	endif()

	# SQL is not needed
	set(DLIB_LINK_WITH_SQLITE3 0)

	# set the root to DLib library
	set(DLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../util/dlib/")

	# include some funny cmake functions from DLib
	include("${DLIB_ROOT}/dlib/cmake/")

	# append to flags definitions
	set(FLAGS_DEF "-USE_DLIB ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_DLIB ${FLAGS_DEF_D}")

	# include header files
	link_directories("${PASCINFERENCE_ROOT}/util/dlib/")
	include_directories("${PASCINFERENCE_ROOT}/util/dlib/")

	# link library
	set(LIBRARIES_DEF "dlib::dlib;${LIBRARIES_DEF}")

	option(DLIB_USE_MKL_FFT OFF)
	option(DLIB_USE_USE_MKL_FFT OFF)
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_DLIB)
	printinfo_onoff("USE_DLIB\t\t\t" "${USE_DLIB}")
	if(${USE_DLIB})
		printinfo(" - DLIB_ROOT\t\t\t" "${DLIB_ROOT}")
	endif()

endmacro()

