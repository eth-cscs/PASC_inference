include_directories("${CMAKE_SOURCE_DIR}/test_util/")

# decide which util to compile
option(TEST_UTIL "TEST_UTIL" OFF)
option(TEST_UTIL_DIFF_NORM_VEC "TEST_UTIL_DIFF_NORM_VEC" OFF)
option(TEST_UTIL_STAT_VEC "TEST_UTIL_STAT_VEC" OFF)
option(TEST_UTIL_PRECISION "TEST_UTIL_PRECISION" OFF)

if(${TEST_UTIL})
	# define shortcut to compile all utils
	getListOfVarsStartingWith("TEST_UTIL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# print info
print("\nUtil tests")
printinfo_onoff(" TEST_UTIL                                 (...)                      " "${TEST_UTIL}")
printinfo_onoff("  TEST_UTIL_DIFF_NORM_VEC                   (ConsoleArg)               " "${TEST_UTIL_DIFF_NORM_VEC}")
printinfo_onoff("  TEST_UTIL_STAT_VEC                        (ConsoleArg)               " "${TEST_UTIL_STAT_VEC}")
printinfo_onoff("  TEST_UTIL_PRECISION                       (ConsoleArg)               " "${TEST_UTIL_PRECISION}")

if(${TEST_UTIL_DIFF_NORM_VEC})
	testadd_executable("test_util/diff_norm_vec.cpp" "diff_norm_vec")
endif()

if(${TEST_UTIL_STAT_VEC})
	testadd_executable("test_util/stat_vec.cpp" "stat_vec")
endif()

if(${TEST_UTIL_PRECISION})
#	if(${USE_CUDA})
#		testadd_executable("test_util/precision.cu" "precision")
#	else()
#		testadd_executable("test_util/precision.cpp" "precision")
#	endif()
endif()

