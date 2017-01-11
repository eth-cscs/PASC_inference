include_directories("${CMAKE_SOURCE_DIR}/test_edf/")

# decide which example to compile
option(TEST_EDF "TEST_EDF" OFF)
option(TEST_EDF2 "TEST_EDF2" OFF)

# print info
printinfo_onoff("TEST_EDF\t\t\t\t" "${TEST_EDF}")
printinfo_onoff("TEST_EDF2\t\t\t\t" "${TEST_EDF2}")

if(${TEST_EDF})
	# this is edf test
	if(${USE_CUDA})
		pascadd_executable("test_edf/test_edf.cu" "test_edf")
	else()
		pascadd_executable("test_edf/test_edf.cpp" "test_edf")
	endif()

	# copy scripts
	make_directory("scripts/test_edf/")
	file(COPY "scripts/" DESTINATION "scripts/test_edf/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_edf/scripts/" DESTINATION "scripts/test_edf/"	FILES_MATCHING PATTERN "*")
			
	# copy data
	file(COPY "test_edf/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")
endif()

if(${TEST_EDF2})
	# this is edf test
	if(${USE_CUDA})
		pascadd_executable("test_edf/test_edf2.cu" "test_edf2")
	else()
		pascadd_executable("test_edf/test_edf2.cpp" "test_edf2")
	endif()

	# copy scripts
	make_directory("scripts/test_edf/")
	file(COPY "scripts/" DESTINATION "scripts/test_edf/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_edf/scripts/" DESTINATION "scripts/test_edf/"	FILES_MATCHING PATTERN "*")
			
	# copy data
	file(COPY "test_edf/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")
endif()
