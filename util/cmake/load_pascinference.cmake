# set variables (options) for whole cmake
option(COMPILE_PASCINFERENCE "COMPILE_PASCINFERENCE" ON)

# include directories to this project
set(PASCINFERENCE_INCLUDE "${PASCINFERENCE_ROOT}/include")
set(PASCINFERENCE_SRC "${PASCINFERENCE_ROOT}/src")
include_directories(${PASCINFERENCE_INCLUDE})
include_directories(${PASCINFERENCE_SRC})

if(${COMPILE_PASCINFERENCE})
	message(STATUS "${Blue}compiling PASCInference library${ColourReset}")

	# from source files create shared library
	if(${USE_GPU})
		CUDA_ADD_LIBRARY(libpascinference ${PASCINFERENCE_SRC}/libpascinference.cu OPTIONS ${LIBRARY_DEFS} SHARED )
	else()
		add_library(libpascinference SHARED "${PASCINFERENCE_SRC}/libpascinference.cpp") # ${LIBRARY_DEFS}" )
	endif()
	
endif()


# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PASCINFERENCE)
	printinfo_green("PASCINFERENCE\t\t\t" "")
	if(${COMPILE_PASCINFERENCE})
		printinfo_green(" - compile\t\t\t" "YES")
	else()
		printinfo_red(" - compile\t\t\t" "NO")
	endif()
	
	printinfo(" - PASCINFERENCE_ROOT\t\t" "${PASCINFERENCE_ROOT}")
	printinfo(" - include\t\t\t" "${PASCINFERENCE_INCLUDE}")
	printinfo(" - src\t\t\t\t" "${PASCINFERENCE_SRC}")

endmacro()

