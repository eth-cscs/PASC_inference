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
set(LIBRARIES_DEF "-lrt;")
set(COMPILE_FIRST "")

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp")

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

link_directories($ENV{LIBRARY_PATH})
include_directories($ENV{CPATH})

macro(PRINT value)
 message("${value}")
endmacro()

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

macro(PRINTINFO_ONOFF name value)
	if(${value})
		message(" ${ColourReset}${name} : ${Green}${value}${Yellow}${ColourReset}")
	else()
		message(" ${ColourReset}${name} : ${Red}${value}${Yellow}${ColourReset}")
	endif()
endmacro()


# add executable file
macro(PASCADD_EXECUTABLE filename outname)
		add_executable(${outname} ${filename})

		# add dependency - build pasc library first and then this exec
		if(${COMPILE_PASCINFERENCE})
			add_dependencies(${outname} libpascinference)
		endif()

		# link pasc inference library
		target_link_libraries(${outname} ${PASCINFERENCE_LIB_LOCATION})

		# link external libraries	
		target_link_libraries(${outname} ${LIBRARIES_DEF})

		# set the name of output file
		set_source_files_properties(${filename}
				COMPILE_FLAGS "${FLAGS_DEF_D}")

		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
#			DEBUG ${CMAKE_CXX_FLAGS_DEBUG}
		)
endmacro()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PASCINFERENCE)
	print("PASCINFERENCE")
	printinfo(" - PASCINFERENCE_ROOT\t\t" "${PASCINFERENCE_ROOT}")
	printinfo(" - include\t\t\t" "${PASCINFERENCE_INCLUDE}")
	printinfo(" - src\t\t\t\t" "${PASCINFERENCE_SRC}")

endmacro()

# get the list of variables which start with some prefix
function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# add library with extension
function (pascinference_add_library _extension _varResult)
    get_cmake_property(_vars VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()


macro(PASCINFERENCE_ADD_LIBRARY extension)
	find_library(PASCINFERENCE_LIB "pascinference_${extension}" "${PASCINFERENCE_ROOT}/build/lib/")
	if(NOT PASCINFERENCE_LIB)
		# if the library doesn't exist, then give error
		message(FATAL_ERROR "\n${Red}PASC_Inference_${extension} not found, did you forget to compile it? : ${PASCINFERENCE_LIB} ${ColourReset}\n")
	else()
		# library found
		message(STATUS "${Blue}PASC_Inference_${extension} library found in: ${PASCINFERENCE_LIB} ${ColourReset}")
	endif()
	
	set(LIBRARIES_DEF "${LIBRARIES_DEF};${PASCINFERENCE_LIB}")
endmacro()
