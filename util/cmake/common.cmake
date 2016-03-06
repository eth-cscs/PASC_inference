# define general functions and variables used in other scripts

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

# here will be loaded options for compilers
set(LIBRARY_DEFS "")
set(DEVICE_DEFS "")

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

