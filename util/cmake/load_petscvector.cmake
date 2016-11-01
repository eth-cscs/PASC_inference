
if(${USE_PETSCVECTOR})
	message(STATUS "${Yellow}loading petscvector library${ColourReset}")
	set(PETSCVECTOR_INCLUDE "${CMAKE_SOURCE_DIR}/../include/algebra/vector/petscvector")
	
	#we have to use Petsc to use PetscVector
	if(NOT ${USE_PETSC})
		message(FATAL_ERROR "${Red}Sorry, you cannot use PetscVector without Petsc! (use -DUSE_PETSC=ON)${ColourReset}")
	endif()
	
	#control existence
	if(NOT EXISTS "${PETSCVECTOR_INCLUDE}/petscvector.h")
		message(FATAL_ERROR "${Red}PetscVector library cannot be found in ${PETSCVECTOR_INCLUDE}${ColourReset}")
	endif()	
	
	#include petscvector
	include_directories(${PETSCVECTOR_INCLUDE})

	# append to flags definitions
	set(FLAGS_DEF "-USE_PETSCVECTOR ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_PETSCVECTOR ${FLAGS_DEF_D}")

endif()


# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PETSCVECTOR)
	printinfo_onoff("USE_PETSCVECTOR\t\t" "${USE_PETSCVECTOR}")
	if(${USE_PETSCVECTOR})
		printinfo(" - include\t\t\t" "${PETSCVECTOR_INCLUDE}")
	endif()
endmacro()

