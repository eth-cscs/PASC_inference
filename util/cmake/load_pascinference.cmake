
# include directories to this project
set(PASCINFERENCE_INCLUDE "${PASCINFERENCE_ROOT}/include")
set(PASCINFERENCE_SRC "${PASCINFERENCE_ROOT}/src")
include_directories(${PASCINFERENCE_INCLUDE})
include_directories(${PASCINFERENCE_SRC})

if(${COMPILE_PASCINFERENCE})
	message(STATUS "${Blue}compiling PASCInference library${ColourReset}")

	# from source files create shared library
	if(${USE_GPU})
		message(STATUS "${Blue}- will be compiled using CUDA compiler${ColourReset}")
		CUDA_ADD_LIBRARY(libpascinference ${PASCINFERENCE_SRC}/libpascinference.cu OPTIONS ${FLAGS_DEF_D} SHARED )
	else()
		message(STATUS "${Blue}- will be compiled using C++ compiler${ColourReset}")
		add_library(libpascinference SHARED "${PASCINFERENCE_SRC}/libpascinference.cpp")
		set_source_files_properties("${PASCINFERENCE_SRC}/libpascinference.cpp"
				COMPILE_FLAGS ${FLAGS_DEF_D})		
	endif()

	set(PASCINFERENCE_LIB_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/liblibpascinference.so")

else()
	# find PascInference library
	find_library(PASCINFERENCE_LIB_LOCATION "libpascinference")
	
endif()

set(LIBRARIES_DEF ${PASCINFERENCE_LIB_LOCATION} "${LIBRARIES_DEF}")

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PASCINFERENCE)
	printinfo_green("PASCINFERENCE\t\t\t" "")
	if(${COMPILE_PASCINFERENCE})
		printinfo_green(" - compile\t\t\t" "YES")
	else()
		printinfo_red(" - compile\t\t\t" "NO")
	endif()
	
	printinfo(" - PASCINFERENCE_ROOT\t\t" "${PASCINFERENCE_ROOT}")
	printinfo(" - library\t\t\t" "${PASCINFERENCE_LIB_LOCATION}")
	printinfo(" - include\t\t\t" "${PASCINFERENCE_INCLUDE}")
	printinfo(" - src\t\t\t\t" "${PASCINFERENCE_SRC}")

endmacro()

