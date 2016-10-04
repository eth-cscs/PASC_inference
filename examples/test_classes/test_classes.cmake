include_directories("${CMAKE_SOURCE_DIR}/test_classes/")

# ----- COMMON -----
option(TEST_COMMON "TEST_COMMON" OFF)
option(TEST_COMMON_CONSOLEARG "TEST_COMMON_CONSOLEARG" OFF)
option(TEST_COMMON_CONSOLEOUTPUT "TEST_COMMON_CONSOLEOUTPUT" OFF)
option(TEST_COMMON_GLOBALMANAGER "TEST_COMMON_GLOBALMANAGER" OFF)
option(TEST_COMMON_LOGGING "TEST_COMMON_LOGGING" OFF)
option(TEST_COMMON_MEMORYCHECK "TEST_COMMON_MEMORYCHECK" OFF)
option(TEST_COMMON_OFFSET "TEST_COMMON_OFFSET" OFF)
option(TEST_COMMON_SHORTINFO "TEST_COMMON_SHORTINFO" OFF)
option(TEST_COMMON_TIMER "TEST_COMMON_TIMER" OFF)
if(${TEST_COMMON})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_COMMON_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# print info
print("\nClasses tests")
printinfo_onoff(" TEST_COMMON                      (...)                " "${TEST_COMMON}")
printinfo_onoff("   TEST_COMMON_CONSOLEARG         (ConsoleArg)         " "${TEST_COMMON_CONSOLEARG}")
printinfo_onoff("   TEST_COMMON_CONSOLEOUTPUT      (ConsoleOutput)      " "${TEST_COMMON_CONSOLEOUTPUT}")
printinfo_onoff("   TEST_COMMON_GLOBALMANAGER      (GlobalManager)      " "${TEST_COMMON_GLOBALMANAGER}")
printinfo_onoff("   TEST_COMMON_LOGGING            (Logging)            " "${TEST_COMMON_LOGGING}")
printinfo_onoff("   TEST_COMMON_MEMORYCHECK        (MemoryCheck)        " "${TEST_COMMON_MEMORYCHECK}")
printinfo_onoff("   TEST_COMMON_OFFSET             (Offset)             " "${TEST_COMMON_OFFSET}")
printinfo_onoff("   TEST_COMMON_SHORTINFO          (Shortinfo)          " "${TEST_COMMON_SHORTINFO}")
printinfo_onoff("   TEST_COMMON_TIMER              (Timer,StackTimer)   " "${TEST_COMMON_TIMER}")


if(${TEST_COMMON_CONSOLEARG})
	# ConsoleArgClass
	pascadd_executable("test_classes/common/test_consolearg.cpp" "test_consolearg")
endif()

if(${TEST_COMMON_CONSOLEOUTPUT})
	# ConsoleOutput
	pascadd_executable("test_classes/common/test_consoleoutput.cpp" "test_consoleoutput")
endif()

if(${TEST_COMMON_GLOBALMANAGER})
	# GlobalManager
	pascadd_executable("test_classes/common/test_globalmanager.cpp" "test_globalmanager")
endif()

if(${TEST_COMMON_LOGGING})
	# Logging
	pascadd_executable("test_classes/common/test_logging.cpp" "test_logging")
endif()

if(${TEST_COMMON_MEMORYCHECK})
	# MemoryCheck
	pascadd_executable("test_classes/common/test_memorycheck.cpp" "test_memorycheck")
endif()

if(${TEST_COMMON_OFFSET})
	# Offset
	pascadd_executable("test_classes/common/test_offset.cpp" "test_offset")
endif()

if(${TEST_COMMON_SHORTINFO})
	# Shortinfo
	pascadd_executable("test_classes/common/test_shortinfo.cpp" "test_shortinfo")
endif()

if(${TEST_COMMON_TIMER})
	# Timer and StackTimer
	pascadd_executable("test_classes/common/test_timer.cpp" "test_timer")
endif()


