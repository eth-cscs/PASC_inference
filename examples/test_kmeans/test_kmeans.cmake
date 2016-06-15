include_directories("${CMAKE_SOURCE_DIR}/test_kmeans/")

# decide which example to compile
option(TEST_KMEANS_LARGE "TEST_KMEANS_LARGE" OFF)

# print info
printinfo_onoff("TEST_KMEANS_LARGE\t\t" "${TEST_KMEANS_LARGE}")

# add example executable files
if(${TEST_KMEANS_LARGE})
	# this is VARX global test
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_kmeans/test_kmeans_large.cu" "test_kmeans_large")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_kmeans/test_kmeans_large.cpp" "test_kmeans_large")
	endif()
endif()

