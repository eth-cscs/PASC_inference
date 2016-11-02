if(${USE_BOOST})
	message(STATUS "${Yellow}loading boost library${ColourReset}")

	find_package( Boost REQUIRED )

	if(Boost_FOUND)
		include_directories( ${Boost_INCLUDE_DIR} )
		link_directories(${Boost_LIBRARY_DIRS})

		# append to flags definitions
		set(FLAGS_DEF "-USE_BOOST ${FLAGS_DEF}")
		set(FLAGS_DEF_D "-DUSE_BOOST ${FLAGS_DEF_D}")

		# todo: here write all used libraries of boost
		set(LIBRARIES_DEF "${LIBRARIES_DEF};boost_system;boost_filesystem;boost_program_options")
	else()
		message(FATAL_ERROR "${Red}Boost library cannot be found.${ColourReset}")
	endif()
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_BOOST)
	printinfo_onoff("USE_BOOST\t\t\t" "${USE_BOOST}")
	if(${USE_BOOST})
		printinfo(" - Boost_INCLUDE_DIR\t\t" "${Boost_INCLUDE_DIR}")
		printinfo(" - Boost_LIBRARY_DIRS\t\t" "${Boost_LIBRARY_DIRS}")
		printinfo(" - Boost_LIBRARIES\t\t" "${Boost_LIBRARIES}")
	endif()
endmacro()

