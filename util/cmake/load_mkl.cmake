# set variables (options) for mkl

if(${USE_MKL})
	message(STATUS "${Blue}loading MKL library${ColourReset}")

	# define variables for mkl include directories
	if(NOT DEFINED ENV{MKLROOT})
		message(FATAL_ERROR "${Red}MKLROOT has to be specified!${ColourReset}")
		return()
	endif()
	set(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
	include_directories(${MKL_INCLUDE_DIR})
	link_directories($ENV{MKLROOT}/lib/intel64)

	# append to flags definitions
	set(FLAGS_DEF "-USE_MKL ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_MKL ${FLAGS_DEF_D}")
	set(LIBRARIES_DEF "mkl_core;mkl_gnu_thread;mkl_rt;${LIBRARIES_DEF}")

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_MKL)
	printinfo_yesno("USE_MKL\t\t\t" "${USE_MKL}")
	if(${USE_MKL})
		printinfo(" - MKL_INCLUDE\t\t" "${MKL_INCLUDE_DIR}")
	endif()

endmacro()

