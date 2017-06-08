include_directories("${CMAKE_SOURCE_DIR}/test_signal/")

# decide which example to compile
option(TEST_SIGNAL "TEST_SIGNAL" OFF)
option(TEST_SIGNAL1D_GENERATE "TEST_SIGNAL1D_GENERATE" OFF)
option(TEST_SIGNAL2D_GENERATE "TEST_SIGNAL2D_GENERATE" OFF)

# print info
print("Signal tests")
printinfo_onoff(" TEST_SIGNAL                                                                          " "${TEST_SIGNAL}")
printinfo_onoff(" TEST_SIGNAL1D_GENERATE                                                               " "${TEST_SIGNAL1D_GENERATE}")
printinfo_onoff(" TEST_SIGNAL2D_GENERATE                                                               " "${TEST_SIGNAL2D_GENERATE}")

if(${TEST_SIGNAL})
	# this is signal processing test
	testadd_executable("test_signal/test_signal.cpp" "test_signal")

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
	file(COPY "test_signal/batch/" DESTINATION "batch/test_signal/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_signal/")
	file(COPY "test_signal/data/" DESTINATION "data/test_signal" FILES_MATCHING PATTERN "*")

endif()

if(${TEST_SIGNAL2D_GENERATE})
	# this is signal processing test
	testadd_executable("test_signal/test_signal2D_generate.cpp" "test_signal2D_generate")

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
