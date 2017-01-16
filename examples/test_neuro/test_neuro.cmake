include_directories("${CMAKE_SOURCE_DIR}/test_neuro/")

# decide which example to compile
option(TEST_NEURO "TEST_NEURO" OFF)

# print info
print("\nNeuro tests")
printinfo_onoff(" TEST_NEURO                                                           " "${TEST_NEURO}")

if(${TEST_NEURO})
	# this is neuro test
	if(${USE_CUDA})
		pascadd_executable("test_neuro/test_neuro.cu" "test_neuro")
	else()
		pascadd_executable("test_neuro/test_neuro.cpp" "test_neuro")
	endif()

	# copy scripts
	make_directory("scripts/test_neuro/")
	file(COPY "scripts/" DESTINATION "scripts/test_neuro/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_neuro/scripts/" DESTINATION "scripts/test_neuro/"	FILES_MATCHING PATTERN "*")
			
	# copy data
	file(COPY "test_neuro/data/" DESTINATION "data" FILES_MATCHING PATTERN "*")
endif()

