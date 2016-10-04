include_directories("${CMAKE_SOURCE_DIR}/test_classes/")

# decide which test to compile
option(TEST_CONSOLEARG "TEST_CONSOLEARG" OFF)
option(TEST_CONSOLEOUTPUT "TEST_CONSOLEOUTPUT" OFF)
option(TEST_GLOBALMANAGER "TEST_GLOBALMANAGER" OFF)
option(TEST_LOGGING "TEST_LOGGING" OFF)
option(TEST_MEMORYCHECK "TEST_MEMORYCHECK" OFF)
option(TEST_OFFSET "TEST_OFFSET" OFF)
option(TEST_SHORTINFO "TEST_SHORTINFO" OFF)
option(TEST_TIMER "TEST_TIMER" OFF)

# print info
print("\nClasses tests")
print(" common")
printinfo_onoff("  TEST_CONSOLEARG         (ConsoleArg)         " "${TEST_CONSOLEARG}")
printinfo_onoff("  TEST_CONSOLEOUTPUT      (ConsoleOutput)      " "${TEST_CONSOLEOUTPUT}")
printinfo_onoff("  TEST_GLOBALMANAGER      (GlobalManager)      " "${TEST_GLOBALMANAGER}")
printinfo_onoff("  TEST_LOGGING            (Logging)            " "${TEST_LOGGING}")
printinfo_onoff("  TEST_MEMORYCHECK        (MemoryCheck)        " "${TEST_MEMORYCHECK}")
printinfo_onoff("  TEST_OFFSET             (Offset)             " "${TEST_OFFSET}")
printinfo_onoff("  TEST_SHORTINFO          (Shortinfo)          " "${TEST_SHORTINFO}")
printinfo_onoff("  TEST_TIMER              (Timer,StackTimer)   " "${TEST_TIMER}")


if(${TEST_CONSOLEARG})
	# ConsoleArgClass
	pascadd_executable("test_classes/test_consolearg.cpp" "test_consolearg")
endif()

if(${TEST_CONSOLEOUTPUT})
	# ConsoleOutput
	pascadd_executable("test_classes/test_consoleoutput.cpp" "test_consoleoutput")
endif()

if(${TEST_GLOBALMANAGER})
	# GlobalManager
	pascadd_executable("test_classes/test_globalmanager.cpp" "test_globalmanager")
endif()

if(${TEST_LOGGING})
	# Logging
	pascadd_executable("test_classes/test_logging.cpp" "test_logging")
endif()

if(${TEST_MEMORYCHECK})
	# MemoryCheck
	pascadd_executable("test_classes/test_memorycheck.cpp" "test_memorycheck")
endif()

if(${TEST_OFFSET})
	# Offset
	pascadd_executable("test_classes/test_offset.cpp" "test_offset")
endif()

if(${TEST_SHORTINFO})
	# Shortinfo
	pascadd_executable("test_classes/test_shortinfo.cpp" "test_shortinfo")
endif()

if(${TEST_TIMER})
	# Timer and StackTimer
	pascadd_executable("test_classes/test_timer.cpp" "test_timer")
endif()


