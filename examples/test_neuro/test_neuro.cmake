include_directories("${CMAKE_SOURCE_DIR}/test_neuro/")

# decide which example to compile
option(TEST_NEURO "TEST_NEURO" OFF)

# print info
print("Neuro tests")
printinfo_onoff(" TEST_NEURO                                                                           " "${TEST_NEURO}")


if(${TEST_NEURO})
	# this is neuro test
	testadd_executable("test_neuro/test_neuro.cpp" "test_neuro")

	# copy scripts
	make_directory("scripts/test_neuro/")
	file(COPY "scripts/" DESTINATION "scripts/test_neuro/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_neuro/scripts/" DESTINATION "scripts/test_neuro/"	FILES_MATCHING PATTERN "*")
			
	# copy data
	make_directory("data/test_neuro/")
	file(COPY "test_neuro/data/" DESTINATION "data/test_neuro/" FILES_MATCHING PATTERN "*")
endif()

