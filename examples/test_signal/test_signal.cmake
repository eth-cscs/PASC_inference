include_directories("${CMAKE_SOURCE_DIR}/test_signal/")

# decide which example to compile
option(TEST_SIGNAL1D "TEST_SIGNAL1D" OFF)
option(TEST_SIGNAL1D_GENERATE "TEST_SIGNAL1D_GENERATE" OFF)

# print info
print("Signal tests")
printinfo_onoff(" TEST_SIGNAL1D                                                                        " "${TEST_SIGNAL1D}")
printinfo_onoff(" TEST_SIGNAL1D_GENERATE                                                               " "${TEST_SIGNAL1D_GENERATE}")

if(${TEST_SIGNAL1D})
	# this is signal processing test
	testadd_executable("test_signal/test_signal1D.cpp" "test_signal1D")

	# copy scripts
	make_directory("scripts/test_signal/")
	file(COPY "scripts/" DESTINATION "scripts/test_signal/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_signal/scripts/" DESTINATION "scripts/test_signal/" FILES_MATCHING PATTERN "*")

	# copy batch scripts
	make_directory("batch/test_signal/")
	file(COPY "test_signal/batch/" DESTINATION "batch/test_signal/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_signal/")
	file(COPY "test_signal/data/" DESTINATION "data/test_signal" FILES_MATCHING PATTERN "*")

endif()

if(${TEST_SIGNAL1D_GENERATE})
	# this is signal processing test
	testadd_executable("test_signal/test_signal1D_generate.cpp" "test_signal1D_generate")

	# copy scripts
	make_directory("scripts/test_signal/")
	file(COPY "scripts/" DESTINATION "scripts/test_signal/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_signal/scripts/" DESTINATION "scripts/test_signal/" FILES_MATCHING PATTERN "*")

	# copy batch scripts
	make_directory("batch/test_signal/")
	file(COPY "test_signal1D/batch/" DESTINATION "batch/test_signal/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_signal/")
	file(COPY "test_signal/data/" DESTINATION "data/test_signal" FILES_MATCHING PATTERN "*")

endif()

