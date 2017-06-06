include_directories("${CMAKE_SOURCE_DIR}/test_image/")

# decide which example to compile
option(TEST_IMAGE "TEST_IMAGE" OFF)

# print info
print("Image tests")
printinfo_onoff(" TEST_IMAGE                                                                           " "${TEST_IMAGE}")


if(${TEST_IMAGE})
	# this is image processing test
	testadd_executable("test_image/test_image.cpp" "test_image")

	# copy scripts
	make_directory("scripts/test_image/")
	file(COPY "scripts/" DESTINATION "scripts/test_image/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_image/scripts/" DESTINATION "scripts/test_image/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_image/")
	file(COPY "test_image/data/" DESTINATION "data/test_image/" FILES_MATCHING PATTERN "*")

endif()

