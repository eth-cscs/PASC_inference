# set variables (options) for whole cmake
option(USE_PETSC "USE_PETSC" ON)

if(${USE_PETSC})
	message(STATUS "${Blue}loading petsc library${ColourReset}")

	# include cmake functions directory
	set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/petsc" ${CMAKE_MODULE_PATH})

	# append to library definitions and device definitions
	set(LIBRARY_DEFS ${LIBRARY_DEFS} -DUSE_PETSC)
	set(DEVICE_DEFS ${DEVICE_DEFS} -USE_PETSC)
	
	# PETSc: include
	if(NOT DEFINED ENV{PETSC_DIR} OR NOT DEFINED ENV{PETSC_INCLUDES})
		message(STATUS "${Blue}PETSC_INCLUDES is not specified, trying to run find_package(PETSc)${ColourReset}")
		find_package(PETSc)
	endif()

	include_directories(${PETSC_INCLUDES})
	
	#include petscvector
#	include_directories(${CMAKE_SOURCE_DIR}/../util/petscvector/include)

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PETSC)
	if(${USE_PETSC})
		printinfo_green("USE_PETSC\t\t\t" "YES")
		printinfo(" - PETSC_DIR\t\t\t" "$ENV{PETSC_DIR}")
		printinfo(" - PETSC_ARCH\t\t\t" "$ENV{PETSC_ARCH}")
		printinfo(" - PETSC_INCLUDES\t\t\t" "$ENV{PETSC_INCLUDES}")
	else()
		printinfo_red("USE_PETSC\t\t" "NO")
	endif()

endmacro()

