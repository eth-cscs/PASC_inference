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

	# include header files
	set(PERMON_INCLUDE_DIR "$ENV{PERMON_DIR}/include")

	if(${USE_CUDA})
		include_directories(${PERMON_INCLUDE_DIR})
		cuda_include_directories(${PERMON_INCLUDE_DIR})
	else()
		include_directories(${PERMON_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_PERMON ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_PERMON ${FLAGS_DEF_D}")

	# include compiled! PERMON library
	link_directories("$ENV{PERMON_DIR}/$ENV{PETSC_ARCH}/lib/")
	find_library(PERMON_LIB_LOCATION "permon" PATH "$ENV{PERMON_DIR}/$ENV{PETSC_ARCH}/lib/")

	#set(PERMON_LIB_LOCATION "${CMAKE_SOURCE_DIR}/../util/metis/build/libmetis/libmetis.so")
	
	if(NOT PERMON_LIB_LOCATION)
		# if the library doesn't exist, then give error
		message(FATAL_ERROR "\n${Red}PERMON library not found, did you forget to compile it? : ${PERMON_LIB_LOCATION} ${ColourReset}\n")
	else()
		# library found
		message(STATUS "${Blue}PERMON library found in: ${PERMON_LIB_LOCATION} ${ColourReset}")
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

