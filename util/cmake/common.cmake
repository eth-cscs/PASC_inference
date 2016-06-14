# define general functions and variables used in other scripts

# include other funny functions
#include(CheckLibraryExists)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

# here will be loaded options for compilers
set(FLAGS_DEF "")
set(FLAGS_DEF_D "")
set(LIBRARIES_DEF "")
set(COMPILE_FIRST "")

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

# define some colors for funny cmake messages
option(CMAKE_USE_COLOR "CMAKE_USE_COLOR" ON)
if(NOT WIN32 AND ${CMAKE_USE_COLOR})
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
endif()

# prepare directory for libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
make_directory(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# print name of variable and value
macro(PRINTINFO name value)
 message(" ${name} : ${Yellow}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_RED name value)
 message(" ${name} : ${Red}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_GREEN name value)
 message(" ${name} : ${Green}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_YESNO name value)
	if(${value})
		message(" ${ColourReset}${name} : ${Green}${value}${Yellow}${ColourReset}")
	else()
		message(" ${ColourReset}${name} : ${Red}${value}${Yellow}${ColourReset}")
	endif()
endmacro()


# add executable file
macro(PASCADD_EXECUTABLE filename outname)
	# choose compiler subject to file extension
	get_filename_component(FILE_EXT ${filename} EXT)

	#remove leading whitespaces
#	string(STRIP ${LIBRARIES_DEF} LIBRARIES_DEF)
	
	if(${FILE_EXT} MATCHES ".cu")
		message(STATUS "${Yellow}adding ${Green}${filename}${Yellow} (CUDA)${ColourReset}")

		# --- compile with CUDA ---
		if(NOT ${USE_CUDA})
			message(FATAL_ERROR "${Red}Cannot compile .cu file without USE_CUDA=ON!${ColourReset}")
		endif()
		
		# add executable file
		cuda_add_executable(${outname} ${filename}
			OPTIONS "${FLAGS_DEF_D} -arch=sm_35 --compiler-options \"${CUDA_CXX_FLAGS}\""
			DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

		# add dependency - build pasc library first and then this exec
		if(${COMPILE_PASCINFERENCE})
			add_dependencies(${outname} libpascinference)
		endif()
	
		# link external libraries	
		target_link_libraries(${outname} ${LIBRARIES_DEF})

		# set the name of output file
		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
		)
		
		# if there are dependencies, then add it
		if(${COMPILE_FIRST})
			add_dependencies(${outname} ${COMPILE_FIRST})
		endif()
	endif()	

	if(${FILE_EXT} MATCHES ".cpp")
		message(STATUS "${Yellow}adding ${Green}${filename}${Yellow} (C++)${ColourReset}")

		# compile with g++

		# add executable file
		add_executable(${outname} ${filename})

		# add dependency - build pasc library first and then this exec
		if(${COMPILE_PASCINFERENCE})
			add_dependencies(${outname} libpascinference)
		endif()

		# link external libraries	
		target_link_libraries(${outname} ${LIBRARIES_DEF})
		

		# set the name of output file
		set_source_files_properties(${filename}
				COMPILE_FLAGS "${FLAGS_DEF_D}")

		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
#			DEBUG ${CMAKE_CXX_FLAGS_DEBUG}
		)
		
		# if there are dependencies, then add it
		if(${COMPILE_FIRST})
			add_dependencies(${outname} ${COMPILE_FIRST})
		endif()
	endif()	
	
	
endmacro()

include(load_cuda) # CUDA
include(load_gpu) # GPU
include(load_petsc) # PETSC
include(load_petscvector) # PetscVector
include(load_boost) # BOOST
include(load_mkl) # MKL
include(load_minlin) # MinLin
include(load_edflib) # EDFlib
include(load_pascinference) # PascInference


