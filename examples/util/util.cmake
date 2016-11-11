include_directories("${CMAKE_SOURCE_DIR}/util/")

# decide which util to compile
option(UTIL "UTIL" OFF)
option(UTIL_DIFF_NORM_VEC "UTIL_DIFF_NORM_VEC" OFF)
option(UTIL_STAT_VEC "UTIL_STAT_VEC" OFF)

if(${UTIL})
	# define shortcut to compile all utils
	getListOfVarsStartingWith("UTIL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# print info
print("\nUtil")
printinfo_onoff(" UTIL                                      (...)                      " "${UTIL}")
printinfo_onoff("  UTIL_DIFF_NORM_VEC                        (ConsoleArg)               " "${UTIL_DIFF_NORM_VEC}")
printinfo_onoff("  UTIL_STAT_VEC                             (ConsoleArg)               " "${UTIL_STAT_VEC}")

if(${UTIL_DIFF_NORM_VEC})
	pascadd_executable("util/diff_norm_vec.cpp" "diff_norm_vec")
endif()

if(${UTIL_STAT_VEC})
	pascadd_executable("util/stat_vec.cpp" "stat_vec")
endif()


