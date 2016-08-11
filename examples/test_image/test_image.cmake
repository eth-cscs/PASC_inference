include_directories("${CMAKE_SOURCE_DIR}/test_image/")

# decide which example to compile
option(TEST_IMAGE "TEST_IMAGE" OFF)
option(TEST_IMAGE_TS "TEST_IMAGE_TS" OFF)

# print info
printinfo_onoff("TEST_IMAGE\t\t\t\t" "${TEST_IMAGE}")
printinfo_onoff("TEST_IMAGE_TS\t\t\t\t" "${TEST_IMAGE_TS}")

if(${TEST_IMAGE})
	# this is image processing test
	if(${USE_CUDA})
		pascadd_executable("test_image/test_image.cu" "test_image")
	else()
		pascadd_executable("test_image/test_image.cpp" "test_image")
	endif()
	
	# copy data
	file(COPY "test_image/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")

	# copy scripts
	make_directory("scripts/test_image/")
	file(COPY "test_image/scripts/" 
		 DESTINATION "scripts/test_image/"
		 FILES_MATCHING PATTERN "*")
endif()

if(${TEST_IMAGE_TS})
	# this is image processing test
	if(${USE_CUDA})
		pascadd_executable("test_image/test_image_ts.cu" "test_image_ts")
	else()
		pascadd_executable("test_image/test_image_ts.cpp" "test_image_ts")
	endif()
	
	# copy data
	file(COPY "test_image/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")
endif()

