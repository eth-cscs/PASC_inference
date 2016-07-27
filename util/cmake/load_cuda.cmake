
if(${USE_CUDA})
	message(STATUS "${Blue}loading CUDA library${ColourReset}")

	include(FindCUDA)
	set(CUDA_PROPAGATE_HOST_FLAGS off) # if flags are passed with -Xcompiler, they also affect NVCC which doesn't understand all g++ flags we use
	set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # without this, cc is used instead of CC and all include paths have to be specified manually
	string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
	set(CUDA_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -lcusparse -Wno-vla ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}") # add flags specific to build type
	string(REPLACE "-std=c++11" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++11 from options
	string(REPLACE "-Wall" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove Wall from options

	set(FLAGS_DEF "-USE_CUDA ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_CUDA ${FLAGS_DEF_D}")
	set(LIBRARIES_DEF "${LIBRARIES_DEF}")
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_CUDA)
	printinfo_onoff("USE_CUDA\t\t\t" "${USE_CUDA}")
	if(${USE_CUDA})
		printinfo(" - CUDA_CXX_FLAGS\t\t" "${CUDA_CXX_FLAGS}")
	endif()

endmacro()

