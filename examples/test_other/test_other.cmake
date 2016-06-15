include_directories("${CMAKE_SOURCE_DIR}/test_other/")

# decide which example to compile
option(TEST_PROJECTION "TEST_PROJECTION" OFF)
option(TEST_GRAPH "TEST_GRAPH" OFF)

# print info
printinfo_onoff("TEST_PROJECTION\t\t" "${TEST_PROJECTION}")
printinfo_onoff("TEST_GRAPH\t\t\t" "${TEST_GRAPH}")

if(${TEST_PROJECTION})
	# this is projection test
	if(${USE_CUDA})
		pascadd_executable("test_other/test_projection.cu" "test_projection")
	else()
		pascadd_executable("test_other/test_projection.cpp" "test_projection")
	endif()
endif()

if(${TEST_GRAPH})
	# this is a test of graphmatrix multiplication
	if(${USE_CUDA})
		pascadd_executable("test_other/test_graph.cu" "test_graph")
	else()
		pascadd_executable("test_other/test_graph.cpp" "test_graph")
	endif()
endif()
