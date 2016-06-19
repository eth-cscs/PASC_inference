include_directories("${CMAKE_SOURCE_DIR}/test_kmeans/")
include_directories("${CMAKE_SOURCE_DIR}/test_kmeans/include")

# decide which example to compile
option(TEST_KMEANS_LARGE "TEST_KMEANS_LARGE" OFF)
option(TEST_KMEANS_LOAD "TEST_KMEANS_LOAD" OFF)
option(TEST_KMEANS_GENERATE "TEST_KMEANS_GENERATE" OFF)

# print info
printinfo_onoff("TEST_KMEANS_LARGE\t\t\t" "${TEST_KMEANS_LARGE}")
printinfo_onoff("TEST_KMEANS_LOAD\t\t\t" "${TEST_KMEANS_LOAD}")
printinfo_onoff("TEST_KMEANS_GENERATE\t\t\t" "${TEST_KMEANS_GENERATE}")

# add example executable files
if(${TEST_KMEANS_LARGE})
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_kmeans/test_kmeans_large.cu" "test_kmeans_large")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_kmeans/test_kmeans_large.cpp" "test_kmeans_large")
	endif()
endif()

if(${TEST_KMEANS_LOAD})
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_kmeans/test_kmeans_load.cu" "test_kmeans_load")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_kmeans/test_kmeans_load.cpp" "test_kmeans_load")
	endif()
	
	# copy data
	file(COPY "test_kmeans/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "data_kmeans_*.bin")
	file(COPY "test_kmeans/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "gamma0_kmeans_*.bin")
		 
endif()

if(${TEST_KMEANS_GENERATE})
	pascadd_executable("test_kmeans/test_kmeans_generate.cpp" "test_kmeans_generate")
endif()
