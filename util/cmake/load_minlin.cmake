# set variables (options) for whole cmake
option(USE_MINLIN "USE_MINLIN" OFF)

if(${USE_MINLIN})
	message(STATUS "${Yellow}loading minlin library${ColourReset}")

	#we have to use Petsc to use PetscVector
	if(NOT ${USE_MKL})
		message(FATAL_ERROR "${Red}Sorry, you cannot use MinLin without MKL!${ColourReset}")
	endif()

	# note: we always need -O3 because minlin_host doesn't compile without it
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

	# MINLIN: define variables for include directories
	set(MINLIN_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../util/minlin/include)
	set(MINLIN_HOST_DEFS
		THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
		THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
		__host__=\ 
		__device__=\ 
		USE_MINLIN)
	set(MINLIN_DEVICE_DEFS # we use -D here because this isn't added automatically
		-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
		-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA
		-DUSE_GPU
		-DUSE_MINLIN)

	include_directories(${MINLIN_INCLUDE})

	if(${USE_GPU})
		cuda_include_directories(${MINLIN_INCLUDE})
	else()
		include_directories(${MINLIN_INCLUDE})
	endif()

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_MINLIN)
	if(${USE_MINLIN})
		printinfo("MINLIN\t\t\t\t" "YES")
		printinfo(" - MINLIN_INCLUDES\t\t" "${MINLIN_INCLUDE}")
	else()
		printinfo_red("MINLIN\t\t\t\t" "NO")
	endif()

endmacro()

