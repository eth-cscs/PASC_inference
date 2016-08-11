include_directories("${CMAKE_SOURCE_DIR}/test_kmeans/")
include_directories("${CMAKE_SOURCE_DIR}/test_kmeans/include")

# decide which example to compile
option(TEST_KMEANS "TEST_KMEANS" OFF)
option(TEST_KMEANS_LOAD "TEST_KMEANS_LOAD" OFF)
option(TEST_KMEANS_GENERATE "TEST_KMEANS_GENERATE" OFF)
option(TEST_KMEANS_SIMPLIFY "TEST_KMEANS_SIMPLIFY" OFF)

# print info
printinfo_onoff("TEST_KMEANS\t\t\t\t" "${TEST_KMEANS}")
printinfo_onoff("TEST_KMEANS_LOAD\t\t\t" "${TEST_KMEANS_LOAD}")
printinfo_onoff("TEST_KMEANS_GENERATE\t\t\t" "${TEST_KMEANS_GENERATE}")
printinfo_onoff("TEST_KMEANS_SIMPLIFY\t\t\t" "${TEST_KMEANS_SIMPLIFY}")

# add example executable files
if(${TEST_KMEANS})
	if(${USE_CUDA})
		# if we know how to compile .cu, then include .cu
		pascadd_executable("test_kmeans/test_kmeans.cu" "test_kmeans")
	else()
		# otherwise compile as .cpp
		pascadd_executable("test_kmeans/test_kmeans.cpp" "test_kmeans")
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
	
	# copy scripts
	make_directory("scripts/test_kmeans/")
	file(COPY "test_kmeans/scripts/" 
		 DESTINATION "scripts/test_kmeans/"
		 FILES_MATCHING PATTERN "*")
	
		 
endif()

# mpiexec -n 1 ./test_kmeans_generate --test_T=10 --test_T=50 --test_T=100 --test_T=500 --test_T=1000 --test_T=5000 --test_T=10000 --test_T=50000 --test_T=100000 --test_T=500000 --test_T=1000000 --test_T=5000000 --test_T=10000000 --test_K=1 --test_K=2 --test_K=3 --test_K=4 --test_K=5
if(${TEST_KMEANS_GENERATE})
	pascadd_executable("test_kmeans/test_kmeans_generate.cpp" "test_kmeans_generate")
endif()

if(${TEST_KMEANS_SIMPLIFY})
	pascadd_executable("test_kmeans/test_kmeans_simplify.cpp" "test_kmeans_simplify")
endif()
