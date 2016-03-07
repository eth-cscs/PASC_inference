# set variables (options) for mkl
option(USE_MKL "USE_MKL" ON)

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
	if(${USE_MKL})
		printinfo_green("USE_MKL\t\t\t" "YES")
		printinfo(" - MKL_INCLUDE\t\t\t" "${MKL_INCLUDE_DIR}")
	else()
		printinfo_red("USE_MKL\t\t\t" "NO")
	endif()

endmacro()

