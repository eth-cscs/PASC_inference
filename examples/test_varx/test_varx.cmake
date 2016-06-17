include_directories("${CMAKE_SOURCE_DIR}/test_varx/")

# decide which example to compile
option(TEST_VARX_LARGE "TEST_VARX_LARGE" OFF)

# print info
printinfo_onoff("TEST_VARX_LARGE\t\t\t" "${TEST_VARX_LARGE}")

# add example executable files
if(${TEST_VARX_LARGE})
	# this is VARX global test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_varx/test_varx_large.cu" "test_varx_large")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_varx/test_varx_large.cpp" "test_varx_large")
	endif()
endif()

