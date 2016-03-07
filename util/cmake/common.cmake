# define general functions and variables used in other scripts

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

# here will be loaded options for compilers
set(FLAGS_DEF "")
set(FLAGS_DEF_D "")
set(LIBRARIES_DEF "")

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

# define some colors for funny cmake messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
endif()

# print name of variable and value
macro(PRINTINFO name value)
 message(" ${name}${Yellow}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_RED name value)
 message(" ${name}${Red}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_GREEN name value)
 message(" ${name}${Green}${value}${ColourReset}")
endmacro()


# add executable file
macro(PASCADD_EXECUTABLE filename outname)
	# choose compiler subject to file extension
	get_filename_component(FILE_EXT ${filename} EXT)
	message(STATUS "${Yellow}adding executable file with extension ${Green}${FILE_EXT}${ColourReset}")

	#remove leading whitespaces
#	string(STRIP ${LIBRARIES_DEF} LIBRARIES_DEF)
	
	if(${FILE_EXT} MATCHES ".cu")
		# --- compile with CUDA ---
		if(NOT ${USE_GPU})
			message(FATAL_ERROR "${Red}Cannot compile .cu file without USE_GPU=ON!${ColourReset}")
		endif()
		
		# add executable file
		cuda_add_executable(${outname} ${filename}
			OPTIONS "${FLAGS_DEF_D} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
			DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	
		# link external libraries	
#		target_link_libraries(${outname} ${LIBRARIES_DEF})
#		target_link_libraries(${outname} "${CMAKE_CURRENT_BINARY_DIR}/pascinference")

		# set the name of output file
		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
		)
	endif()	

	if(${FILE_EXT} MATCHES ".cpp")
		# compile with g++

		# add executable file
		add_executable(${outname} ${filename})

		# link external libraries	
		target_link_libraries(${outname} ${LIBRARIES_DEF})

		# set the name of output file
		set_source_files_properties(${filename}
				COMPILE_FLAGS "${FLAGS_DEF_D}")

		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
#			DEBUG ${CMAKE_CXX_FLAGS_DEBUG}
		)

	endif()	
	
	
endmacro()


