include_directories("${CMAKE_SOURCE_DIR}/util/")

# decide which util to compile
option(UTIL "UTIL" OFF)
option(UTIL_DIFF_NORM_VEC "UTIL_DIFF_NORM_VEC" OFF)
option(UTIL_STAT_VEC "UTIL_STAT_VEC" OFF)
option(UTIL_PRINT_VEC "UTIL_PRINT_VEC" OFF)
option(UTIL_PROCESS_LOG "UTIL_PROCESS_LOG" OFF)

if(${UTIL})
	# define shortcut to compile all utils
	getListOfVarsStartingWith("UTIL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# print info
print("Utils")
printinfo_onoff(" UTIL                                                  (...)                          " "${UTIL}")
printinfo_onoff("   UTIL_DIFF_NORM_VEC                                    (ConsoleArg)                 " "${UTIL_DIFF_NORM_VEC}")
printinfo_onoff("   UTIL_STAT_VEC                                         (ConsoleArg)                 " "${UTIL_STAT_VEC}")
printinfo_onoff("   UTIL_PRINT_VEC                                        (ConsoleArg)                 " "${UTIL_PRINT_VEC}")
printinfo_onoff("   UTIL_PROCESS_LOG                                      (ConsoleArg)                 " "${UTIL_PROCESS_LOG}")

if(${UTIL_DIFF_NORM_VEC})
	testadd_executable("util/util_diff_norm_vec.cpp" "util_diff_norm_vec")
endif()

if(${UTIL_STAT_VEC})
	testadd_executable("util/util_stat_vec.cpp" "util_stat_vec")
endif()

if(${UTIL_PRINT_VEC})
	testadd_executable("util/util_print_vec.cpp" "util_print_vec")
endif()

if(${UTIL_PROCESS_LOG})
	testadd_executable("util/util_process_log.cpp" "util_process_log")
endif()


