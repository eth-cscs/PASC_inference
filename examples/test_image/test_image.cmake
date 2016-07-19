include_directories("${CMAKE_SOURCE_DIR}/test_image/")

# decide which example to compile
option(TEST_IMAGE "TEST_IMAGE" OFF)

# print info
printinfo_onoff("TEST_IMAGE\t\t\t\t" "${TEST_IMAGE}")

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
endif()

