# set variables (options) for minlin
option(USE_MINLIN "USE_MINLIN" ON)

if(${USE_MINLIN})
	message(STATUS "${Yellow}loading minlin library${ColourReset}")

	#we have to use MKL to use MinLin
	if(NOT ${USE_MKL})
		message(FATAL_ERROR "${Red}Sorry, you cannot use MinLin without MKL! (use -DUSE_MKL=ON)${ColourReset}")
	endif()

	#we have to use boost to use MinLin
	if({USE_MKL})
		message(FATAL_ERROR "${Red}Sorry, you cannot use MinLin without BOOST! (use -DUSE_BOOST=ON)${ColourReset}")
	endif()

	# note: we always need -O3 because minlin_host doesn't compile without it
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

	# MINLIN: define variables for include directories
	set(MINLIN_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../util/minlin/include)

	#control existence
	if(NOT EXISTS "${MINLIN_INCLUDE}/minlin/minlin.h")
		message(FATAL_ERROR "${Red}MinLin library cannot be found in ${MINLIN_INCLUDE}${ColourReset}")
	endif()	

	if(${USE_CUDA})
		include_directories(${MINLIN_INCLUDE})
		cuda_include_directories(${MINLIN_INCLUDE})
	else()
		include_directories(${MINLIN_INCLUDE})
	endif()

	# append to flags definitions
	set(FLAGS_DEF "-USE_MINLIN ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_MINLIN ${FLAGS_DEF_D}")

	if(${USE_GPU})
		set(FLAGS_DEF "THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA ${FLAGS_DEF}")
		set(FLAGS_DEF_D "-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA ${FLAGS_DEF_D}")
	else()
		set(FLAGS_DEF "THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP __host__= __device__= ${FLAGS_DEF}")
		set(FLAGS_DEF_D "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP -D__host__= -D__device__= ${FLAGS_DEF_D}")
	endif()

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_MINLIN)
	if(${USE_MINLIN})
		printinfo_green("MINLIN\t\t\t\t" "YES")
		printinfo(" - MINLIN_INCLUDES\t\t" "${MINLIN_INCLUDE}")
	else()
		printinfo_red("MINLIN\t\t\t\t" "NO")
	endif()

endmacro()

