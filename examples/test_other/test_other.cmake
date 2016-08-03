include_directories("${CMAKE_SOURCE_DIR}/test_other/")

# decide which example to compile
option(TEST_PROJECTION "TEST_PROJECTION" OFF)
option(TEST_METIS "TEST_METIS" OFF)
option(TEST_DECOMPOSITION "TEST_DECOMPOSITION" OFF)
option(TEST_MAT_TO_DENSE "TEST_MAT_TO_DENSE" OFF)
option(TEST_MAT_SEQ_VS_MPI "TEST_MAT_SEQ_VS_MPI" OFF)

# print info
printinfo_onoff("TEST_PROJECTION\t\t\t" "${TEST_PROJECTION}")
printinfo_onoff("TEST_METIS\t\t\t\t" "${TEST_METIS}")
printinfo_onoff("TEST_DECOMPOSITION\t\t\t" "${TEST_DECOMPOSITION}")
printinfo_onoff("TEST_MAT_TO_DENSE\t\t\t" "${TEST_MAT_TO_DENSE}")
printinfo_onoff("TEST_MAT_SEQ_VS_MPI\t\t\t" "${TEST_MAT_SEQ_VS_MPI}")

if(${TEST_PROJECTION})
	# this is projection test
	if(${USE_CUDA})
		pascadd_executable("test_other/test_projection.cu" "test_projection")
	else()
		pascadd_executable("test_other/test_projection.cpp" "test_projection")
	endif()
endif()

if(${TEST_METIS})
	# this is test of graph decomposition
	if(${USE_CUDA})
		pascadd_executable("test_other/test_metis.cu" "test_metis")
	else()
		pascadd_executable("test_other/test_metis.cpp" "test_metis")
	endif()

	# copy data
	file(COPY "test_other/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")
endif()

if(${TEST_DECOMPOSITION})
	# this is test of decomposition
	if(${USE_CUDA})
		pascadd_executable("test_other/test_decomposition.cu" "test_decomposition")
	else()
		pascadd_executable("test_other/test_decomposition.cpp" "test_decomposition")
	endif()

	# copy data
	file(COPY "test_other/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")
endif()

if(${TEST_MAT_TO_DENSE})
	# this is a test of matrix multiplication
	if(${USE_CUDA})
		pascadd_executable("test_other/test_mat_to_dense.cu" "test_mat_to_dense")
	else()
		pascadd_executable("test_other/test_mat_to_dense.cpp" "test_mat_to_dense")
	endif()

	# copy data
	file(COPY "test_other/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")
endif()

if(${TEST_MAT_SEQ_VS_MPI})
	if(${USE_CUDA})
		pascadd_executable("test_other/test_mat_seq_vs_mpi.cu" "test_mat_seq_vs_mpi")
	else()
		pascadd_executable("test_other/test_mat_seq_vs_mpi.cpp" "test_mat_seq_vs_mpi")
	endif()

	# copy data
	file(COPY "test_other/data/" 
		 DESTINATION "data"
		 FILES_MATCHING PATTERN "*")	
endif()
