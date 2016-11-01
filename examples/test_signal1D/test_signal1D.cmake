include_directories("${CMAKE_SOURCE_DIR}/test_signal1D/")

# decide which example to compile
option(TEST_SIGNAL1D "TEST_SIGNAL1D" OFF)

# print info
print("\nSignal1D tests")
printinfo_onoff(" TEST_SIGNAL1D                                                        " "${TEST_SIGNAL1D}")

if(${TEST_SIGNAL1D})
	# this is signal processing test
	if(${USE_CUDA})
		pascadd_executable("test_signal1D/test_signal1D.cu" "test_signal1D")
	else()
		pascadd_executable("test_signal1D/test_signal1D.cpp" "test_signal1D")
	endif()

	# copy scripts
	make_directory("scripts/test_signal1D/")
	file(COPY "scripts/" DESTINATION "scripts/test_signal1D/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_signal1D/scripts/" DESTINATION "scripts/test_signal1D/" FILES_MATCHING PATTERN "*")
	
	# copy data
	file(COPY "test_signal1D/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")

endif()

