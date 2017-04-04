
if(${USE_CUDA})
	message(STATUS "${Blue}loading CUDA library${ColourReset}")

	include(FindCUDA)
	set(CUDA_PROPAGATE_HOST_FLAGS off) # if flags are passed with -Xcompiler, they also affect NVCC which doesn't understand all g++ flags we use
	set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # without this, cc is used instead of CC and all include paths have to be specified manually
	string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
	set(CUDA_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lssl -lcusparse -Wno-vla ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}") # add flags specific to build type
	string(REPLACE "-std=c++11" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++11 from options
	string(REPLACE "-std=c++0x" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove C++0x from options
	string(REPLACE "-Wall" "" CUDA_CXX_FLAGS ${CUDA_CXX_FLAGS}) # remove Wall from options

	set(FLAGS_DEF "-USE_CUDA ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_CUDA ${FLAGS_DEF_D}")
	set(LIBRARIES_DEF "ssl;cusparse;${LIBRARIES_DEF}")
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_CUDA)
	printinfo_onoff("USE_CUDA\t\t\t" "${USE_CUDA}")
	if(${USE_CUDA})
		printinfo(" - CUDA_CXX_FLAGS\t\t" "${CUDA_CXX_FLAGS}")
	endif()

endmacro()

macro(PASCADD_EXECUTABLE_CUDA filename outname)
	# --- compile with CUDA ---
	if(NOT ${USE_CUDA})
		message(FATAL_ERROR "${Red}Cannot compile .cu file without USE_CUDA=ON!${ColourReset}")
	endif()
		
	# add executable file
	cuda_add_executable(${outname} ${filename}
			OPTIONS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\""
			DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

	# link external libraries	
	target_link_libraries(${outname} ${LIBRARIES_DEF})

	# set the name of output file
	set_target_properties(${outname} PROPERTIES
		OUTPUT_NAME ${outname}
	)
endmacro()
