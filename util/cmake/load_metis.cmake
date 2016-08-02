# http://glaros.dtc.umn.edu/gkhome/metis/metis/download
# 
# compile using: 
# cd util/metis/build
# cmake -DSHARED=ON ..
# make
#

if(${USE_METIS})
	message(STATUS "${Yellow}loading METIS${ColourReset}")

	set(METIS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../util/metis/include")

	if(${USE_CUDA})
		include_directories(${METIS_INCLUDE_DIR})
		cuda_include_directories(${METIS_INCLUDE_DIR})
	else()
		include_directories(${METIS_INCLUDE_DIR})
	endif()

	set(FLAGS_DEF "-USE_METIS ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_METIS ${FLAGS_DEF_D}")

	# include compiled! METIS library
	link_directories("${CMAKE_SOURCE_DIR}/../util/metis/build/libmetis/")
	find_library(METIS_LIB_LOCATION "metis" PATH "${CMAKE_SOURCE_DIR}/../util/metis/build/libmetis/")

	#set(METIS_LIB_LOCATION "${CMAKE_SOURCE_DIR}/../util/metis/build/libmetis/libmetis.so")
	
	if(NOT METIS_LIB_LOCATION)
		# if the library doesn't exist, then give error
		message(FATAL_ERROR "\n${Red}METIS library not found, did you forget to compile it? : ${METIS_LIB_LOCATION} ${ColourReset}\n")
	else()
		# library found
		message(STATUS "${Blue}METIS library found in: ${METIS_LIB_LOCATION} ${ColourReset}")

	endif()

	# link library
	set(LIBRARIES_DEF "metis;${LIBRARIES_DEF}")

endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_METIS)
	printinfo_onoff("USE_METIS\t\t\t" "${USE_METIS}")
	if(${USE_METIS})
		printinfo(" - METIS_INCLUDE_DIR\t\t" "${METIS_INCLUDE_DIR}")
		printinfo(" - METIS_LIB_LOCATION\t\t" "${METIS_LIB_LOCATION}")
	endif()
endmacro()

