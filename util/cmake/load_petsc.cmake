
if(${USE_PETSC})
	message(STATUS "${Blue}loading petsc library${ColourReset}")

	# include cmake functions directory
	set(CMAKE_MODULE_PATH "${PASCINFERENCE_ROOT}/util/cmake/petsc" ${CMAKE_MODULE_PATH})
	
	# PETSc: include
	if(${FIND_PETSC})
		message(STATUS "${Blue}trying to run find_package(PETSc)${ColourReset}")
		# magic function from Jed Brown
		find_package(PETSc)

		include_directories(${PETSC_INCLUDES})

#		set(CMAKE_CXX_COMPILER "${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx")
#		set(CMAKE_CXX_COMPILER "mpicxx")

	endif()

	# append to flags definitions
	set(FLAGS_DEF "-USE_PETSC ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_PETSC ${FLAGS_DEF_D}")
	set(LIBRARIES_DEF ${PETSC_LIBRARIES} ${LIBRARIES_DEF})

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PETSC)
	if(${USE_PETSC})
		printinfo_green("USE_PETSC\t\t\t" "YES")
		printinfo(" - PETSC_DIR\t\t\t" "$ENV{PETSC_DIR}")
		printinfo(" - PETSC_ARCH\t\t\t" "$ENV{PETSC_ARCH}")
		printinfo(" - PETSC_INCLUDES\t\t" "${PETSC_INCLUDES}")
		printinfo(" - PETSC_LIBRARIES\t\t" "${PETSC_LIBRARIES}")
		printinfo(" - FIND_PETSC\t\t\t" "${FIND_PETSC}")
	else()
		printinfo_red("USE_PETSC\t\t\t" "NO")
	endif()

endmacro()

