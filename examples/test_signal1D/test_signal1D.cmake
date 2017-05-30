include_directories("${CMAKE_SOURCE_DIR}/test_signal1D/")

# decide which example to compile
option(TEST_SIGNAL1D "TEST_SIGNAL1D" OFF)
option(TEST_SIGNAL1D_GENERATE "TEST_SIGNAL1D_GENERATE" OFF)

# print info
print("Signal1D tests")
printinfo_onoff(" TEST_SIGNAL1D                                                                        " "${TEST_SIGNAL1D}")
printinfo_onoff(" TEST_SIGNAL1D_GENERATE                                                               " "${TEST_SIGNAL1D_GENERATE}")

if(${TEST_SIGNAL1D})
	# this is signal processing test
	testadd_executable("test_signal1D/test_signal1D.cpp" "test_signal1D")

	# copy scripts
	make_directory("scripts/test_signal1D/")
	file(COPY "scripts/" DESTINATION "scripts/test_signal1D/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_signal1D/scripts/" DESTINATION "scripts/test_signal1D/" FILES_MATCHING PATTERN "*")
	
	# copy data
	file(COPY "test_signal1D/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")

endif()

if(${TEST_SIGNAL1D_GENERATE})
	# this is signal processing test
	testadd_executable("test_signal1D/test_signal1D_generate.cpp" "test_signal1D_generate")

	# copy scripts
	make_directory("scripts/test_signal1D/")
	file(COPY "scripts/" DESTINATION "scripts/test_signal1D/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_signal1D/scripts/" DESTINATION "scripts/test_signal1D/" FILES_MATCHING PATTERN "*")
	
	# copy data
	file(COPY "test_signal1D/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")

endif()

