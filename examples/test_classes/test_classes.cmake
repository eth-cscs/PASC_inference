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

# ----- ALGEBRA -----
option(TEST_ALGEBRA "TEST_ALGEBRA" OFF)
option(TEST_ALGEBRA_BGMGRAPH "TEST_ALGEBRA_BGMGRAPH" OFF)
option(TEST_ALGEBRA_BGMGRAPHGRID1D "TEST_ALGEBRA_BGMGRAPHGRID1D" OFF)
option(TEST_ALGEBRA_BGMGRAPHGRID2D "TEST_ALGEBRA_BGMGRAPHGRID2D" OFF)
option(TEST_ALGEBRA_BLOCKDIAGMATRIX "TEST_ALGEBRA_BLOCKDIAGMATRIX" OFF)
option(TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX "TEST_ALGEBRA_BLOCKSPARSEMATRIX" OFF)
option(TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX "TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX" OFF)
option(TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX "TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX" OFF)
option(TEST_ALGEBRA_DECOMPOSITION "TEST_ALGEBRA_DECOMPOSITION" OFF)
option(TEST_ALGEBRA_FILECRSMATRIX "TEST_ALGEBRA_FILECRSMATRIX" OFF)
option(TEST_ALGEBRA_PETSCVECTOR "TEST_ALGEBRA_PETSCVECTOR" OFF)
option(TEST_ALGEBRA_LOCALDENSEMATRIX "TEST_ALGEBRA_LOCALDENSEMATRIX" OFF)
option(TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL "TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL" OFF)
if(${TEST_ALGEBRA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_ALGEBRA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# print info
print("\nClasses tests")
printinfo_onoff(" TEST_COMMON                               (...)                      " "${TEST_COMMON}")
printinfo_onoff("   TEST_COMMON_CONSOLEARG                  (ConsoleArg)                " "${TEST_COMMON_CONSOLEARG}")
printinfo_onoff("   TEST_COMMON_CONSOLEOUTPUT               (ConsoleOutput)             " "${TEST_COMMON_CONSOLEOUTPUT}")
printinfo_onoff("   TEST_COMMON_GLOBALMANAGER               (GlobalManager)             " "${TEST_COMMON_GLOBALMANAGER}")
printinfo_onoff("   TEST_COMMON_LOGGING                     (Logging)                   " "${TEST_COMMON_LOGGING}")
printinfo_onoff("   TEST_COMMON_MEMORYCHECK                 (MemoryCheck)               " "${TEST_COMMON_MEMORYCHECK}")
printinfo_onoff("   TEST_COMMON_OFFSET                      (Offset)                    " "${TEST_COMMON_OFFSET}")
printinfo_onoff("   TEST_COMMON_SHORTINFO                   (Shortinfo)                 " "${TEST_COMMON_SHORTINFO}")
printinfo_onoff("   TEST_COMMON_TIMER                       (Timer,StackTimer)          " "${TEST_COMMON_TIMER}")
printinfo_onoff(" TEST_ALGEBRA                              (...)                      " "${TEST_ALGEBRA}")
printinfo_onoff("   TEST_ALGEBRA_BGMGRAPH                   (BGMGraph)                  " "${TEST_ALGEBRA_BGMGRAPH}")
printinfo_onoff("   TEST_ALGEBRA_BGMGRAPHGRID1D             (BGMGraphGrid1D)            " "${TEST_ALGEBRA_BGMGRAPHGRID1D}")
printinfo_onoff("   TEST_ALGEBRA_BGMGRAPHGRID2D             (BGMGraphGrid2D)            " "${TEST_ALGEBRA_BGMGRAPHGRID2D}")
printinfo_onoff("   TEST_ALGEBRA_BLOCKDIAGMATRIX            (BlockDiagMatrix)           " "${TEST_ALGEBRA_BLOCKDIAGMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX     (BlockGraphSparseMatrix)    " "${TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX     (BlockLaplaceFreeMatrix)    " "${TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX   (BlockLaplaceSparseMatrix)  " "${TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_DECOMPOSITION              (Decomposition)             " "${TEST_ALGEBRA_DECOMPOSITION}")
printinfo_onoff("   TEST_ALGEBRA_FILECRSMATRIX              (FileCRSMatrix)             " "${TEST_ALGEBRA_FILECRSMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_PETSCVECTOR                (PetscVector)               " "${TEST_ALGEBRA_PETSCVECTOR}")
printinfo_onoff("   TEST_ALGEBRA_LOCALDENSEMATRIX           (LocalDenseMatrix)          " "${TEST_ALGEBRA_LOCALDENSEMATRIX}")
printinfo_onoff("   TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL    (SimplexFeasibleSet_Local)  " "${TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL}")

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

if(${TEST_ALGEBRA_BGMGRAPH})
	# BGMGraph
	if(${USE_CUDA})
		pascadd_executable("test_classes/algebra/test_bgmgraph.cu" "test_bgmgraph")
	else()
		pascadd_executable("test_classes/algebra/test_bgmgraph.cpp" "test_bgmgraph")
	endif()
	
	# copy data with sample graphs
	file(COPY "test_classes/data/test_algebra_bgmgraph/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

if(${TEST_ALGEBRA_BGMGRAPHGRID1D})
	# BGMGraphGrid1D
	if(${USE_CUDA})
		pascadd_executable("test_classes/algebra/test_bgmgraphgrid1D.cu" "test_bgmgraphgrid1D")
	else()
		pascadd_executable("test_classes/algebra/test_bgmgraphgrid1D.cpp" "test_bgmgraphgrid1D")
	endif()
endif()

if(${TEST_ALGEBRA_BGMGRAPHGRID2D})
	# BGMGraphGrid2D
	if(${USE_CUDA})
		pascadd_executable("test_classes/algebra/test_bgmgraphgrid2D.cu" "test_bgmgraphgrid2D")
	else()
		pascadd_executable("test_classes/algebra/test_bgmgraphgrid2D.cpp" "test_bgmgraphgrid2D")
	endif()
endif()

