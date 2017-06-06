include_directories("${CMAKE_SOURCE_DIR}/test_entropy/")

# decide which example to compile
option(TEST_ENTROPY "TEST_ENTROPY" OFF)

# print info
print("Entropy tests")
printinfo_onoff(" TEST_ENTROPY                                                                         " "${TEST_ENTROPY}")

if(${TEST_ENTROPY})
	# this is signal processing test
	testadd_executable("test_entropy/test_entropy.cpp" "test_entropy")

	# copy scripts
	make_directory("scripts/test_entropy/")
	file(COPY "scripts/" DESTINATION "scripts/test_entropy/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_entropy/scripts/" DESTINATION "scripts/test_entropy/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_entropy/")
	file(COPY "test_entropy/data/" DESTINATION "data/test_entropy/" FILES_MATCHING PATTERN "*")

endif()

