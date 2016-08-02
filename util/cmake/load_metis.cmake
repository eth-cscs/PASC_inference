if(${USE_METIS})
	message(STATUS "${Yellow}loading METIS${ColourReset}")

	set(METIS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../util/metis/include)

	if(${USE_CUDA})
		include_directories(${METIS_INCLUDE_DIR})
		cuda_include_directories(${METIS_INCLUDE_DIR})
	else()
		include_directories(${METIS_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_METIS ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_METIS ${FLAGS_DEF_D}")

endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_METIS)
	printinfo_onoff("USE_METIS\t\t\t" "${USE_METIS}")
	if(${USE_EDFLIB})
		printinfo(" - METIS_INCLUDE_DIR\t\t" "${METIS_INCLUDE_DIR}")
	endif()
endmacro()

