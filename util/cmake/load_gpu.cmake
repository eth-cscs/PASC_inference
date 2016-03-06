# set variables (options) for whole cmake
option(USE_GPU "USE_GPU" OFF)

if(${USE_GPU})
	message(STATUS "${Blue}loading CUDA library${ColourReset}")

	include(FindCUDA)
	set(CUDA_PROPAGATE_HOST_FLAGS off) # if flags are passed with -Xcompiler, they also affect NVCC which doesn't understand all g++ flags we use
	set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # without this, cc is used instead of CC and all include paths have to be specified manually
	string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
	set(CUDA_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-vla ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}") # add flags specific to build type
	string(REPLACE "-std=c++11" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++11 from options

	set(LIBRARY_DEFS ${LIBRARY_DEFS} -DUSE_GPU)
	set(DEVICE_DEFS ${DEVICE_DEFS} -USE_GPU)

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_GPU)
	if(${USE_GPU})
		printinfo_green("USE_GPU\t\t\t" "YES")
		printinfo(" - CUDA_CXX_FLAGS\t\t" "${CUDA_CXX_FLAGS}")
	else()
		printinfo_red("USE_GPU\t\t\t" "NO")
	endif()

endmacro()

