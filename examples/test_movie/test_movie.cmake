include_directories("${CMAKE_SOURCE_DIR}/test_movie/")

# decide which example to compile
option(TEST_MOVIE "TEST_MOVIE" OFF)
option(TEST_MOVIE_GENERATE_DOT "TEST_MOVIE_GENERATE_DOT" OFF)
option(TEST_MOVIE_GENERATE_3DOT "TEST_MOVIE_GENERATE_3DOT" OFF)

# print info
print("Movie tests")
printinfo_onoff(" TEST_MOVIE                                                                           " "${TEST_MOVIE}")
printinfo_onoff(" TEST_MOVIE_GENERATE_DOT                                                              " "${TEST_MOVIE_GENERATE_DOT}")
printinfo_onoff(" TEST_MOVIE_GENERATE_3DOT                                                             " "${TEST_MOVIE_GENERATE_3DOT}")


if(${TEST_MOVIE})
	# this is movie processing test
	testadd_executable("test_movie/test_movie.cpp" "test_movie")

	# copy scripts
	make_directory("scripts/test_movie/")
	file(COPY "scripts/" DESTINATION "scripts/test_movie/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_movie/scripts/" DESTINATION "scripts/test_movie/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_movie/")
	file(COPY "test_movie/data/" DESTINATION "data/test_movie/" FILES_MATCHING PATTERN "*")

endif()

if(${TEST_MOVIE_GENERATE_DOT})
	# this is movie processing test
	testadd_executable("test_movie/test_movie_generate_dot.cpp" "test_movie_generate_dot")

	# copy scripts
	make_directory("scripts/test_movie/")
	file(COPY "scripts/" DESTINATION "scripts/test_movie/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_movie/scripts/" DESTINATION "scripts/test_movie/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_movie/")
	file(COPY "test_movie/data/" DESTINATION "data/test_movie/" FILES_MATCHING PATTERN "*")

endif()

if(${TEST_MOVIE_GENERATE_3DOT})
	# this is colour movie processing test
	testadd_executable("test_movie/test_movie_generate_3dot.cpp" "test_movie_generate_3dot")

	# copy scripts
	make_directory("scripts/test_movie/")
	file(COPY "scripts/" DESTINATION "scripts/test_movie/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_movie/scripts/" DESTINATION "scripts/test_movie/" FILES_MATCHING PATTERN "*")
	
	# copy data
	make_directory("data/test_movie/")
	file(COPY "test_movie/data/" DESTINATION "data/test_movie/" FILES_MATCHING PATTERN "*")

endif()
