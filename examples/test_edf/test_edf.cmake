include_directories("${CMAKE_SOURCE_DIR}/test_edf/")

# decide which example to compile
option(TEST_EDF "TEST_EDF" OFF)

# print info
printinfo_onoff("TEST_EDF\t\t\t" "${TEST_EDF}")

if(${TEST_EDF})
	# this is edf test
	if(${USE_CUDA})
		pascadd_executable("test_edf.cu" "test_edf")
	else()
		pascadd_executable("test_edf.cpp" "test_edf")
	endif()
endif()

