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
option(TEST_ALGEBRA_DOT "TEST_ALGEBRA_DOT" OFF)
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

# ----- DATA -----
option(TEST_DATA "TEST_DATA" OFF)
option(TEST_DATA_DIAG "TEST_DATA_DIAG" OFF)
option(TEST_DATA_EDF "TEST_DATA_EDF" OFF)
option(TEST_DATA_IMAGE "TEST_DATA_IMAGE" OFF)
option(TEST_DATA_KMEANS "TEST_DATA_KMEANS" OFF)
option(TEST_DATA_QP "TEST_DATA_QP" OFF)
option(TEST_DATA_SIGNAL1D "TEST_DATA_SIGNAL1D" OFF)
option(TEST_DATA_SIMPLE "TEST_DATA_SIMPLE" OFF)
option(TEST_DATA_TS "TEST_DATA_TS" OFF)
if(${TEST_DATA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_DATA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- MODEL -----
option(TEST_MODEL "TEST_MODEL" OFF)
option(TEST_MODEL_GRAPHH1FEM "TEST_MODEL_GRAPHH1FEM" OFF)
option(TEST_MODEL_KMEANSH1FEM "TEST_MODEL_KMEANSH1FEM" OFF)
if(${TEST_MODEL})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_MODEL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- SOLVER -----
option(TEST_SOLVER "TEST_SOLVER" OFF)
option(TEST_SOLVER_CGQP "TEST_SOLVER_CGQP" OFF)
option(TEST_SOLVER_DIAG "TEST_SOLVER_DIAG" OFF)
option(TEST_SOLVER_MULTICG "TEST_SOLVER_MULTICG" OFF)
option(TEST_SOLVER_SIMPLE "TEST_SOLVER_SIMPLE" OFF)
option(TEST_SOLVER_SPGQP "TEST_SOLVER_SPGQP" OFF)
if(${TEST_SOLVER})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SOLVER_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- DLIB -----
option(TEST_DLIB "TEST_DLIB" OFF)
option(TEST_DLIB_ANNA "TEST_DLIB_ANNA" OFF)
option(TEST_DLIB_INTEGRAL "TEST_DLIB_INTEGRAL" OFF)
option(TEST_DLIB_GUI "TEST_DLIB_GUI" OFF)
if(${TEST_DLIB})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_DLIB_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# print info
print("\nClasses tests")
printinfo_onoff(" TEST_COMMON                               (...)                      " "${TEST_COMMON}")
printinfo_onoff("  TEST_COMMON_CONSOLEARG                    (ConsoleArg)               " "${TEST_COMMON_CONSOLEARG}")
printinfo_onoff("  TEST_COMMON_CONSOLEOUTPUT                 (ConsoleOutput)            " "${TEST_COMMON_CONSOLEOUTPUT}")
printinfo_onoff("  TEST_COMMON_GLOBALMANAGER                 (GlobalManager)            " "${TEST_COMMON_GLOBALMANAGER}")
printinfo_onoff("  TEST_COMMON_LOGGING                       (Logging)                  " "${TEST_COMMON_LOGGING}")
printinfo_onoff("  TEST_COMMON_MEMORYCHECK                   (MemoryCheck)              " "${TEST_COMMON_MEMORYCHECK}")
printinfo_onoff("  TEST_COMMON_OFFSET                        (Offset)                   " "${TEST_COMMON_OFFSET}")
printinfo_onoff("  TEST_COMMON_SHORTINFO                     (Shortinfo)                " "${TEST_COMMON_SHORTINFO}")
printinfo_onoff("  TEST_COMMON_TIMER                         (Timer,StackTimer)         " "${TEST_COMMON_TIMER}")
printinfo_onoff(" TEST_ALGEBRA                              (...)                      " "${TEST_ALGEBRA}")
printinfo_onoff("  TEST_ALGEBRA_DOT                          (dot product)              " "${TEST_ALGEBRA_DOT}")
printinfo_onoff("  TEST_ALGEBRA_BGMGRAPH                     (BGMGraph)                 " "${TEST_ALGEBRA_BGMGRAPH}")
printinfo_onoff("  TEST_ALGEBRA_BGMGRAPHGRID1D               (BGMGraphGrid1D)           " "${TEST_ALGEBRA_BGMGRAPHGRID1D}")
printinfo_onoff("  TEST_ALGEBRA_BGMGRAPHGRID2D               (BGMGraphGrid2D)           " "${TEST_ALGEBRA_BGMGRAPHGRID2D}")
#printinfo_onoff("  TEST_ALGEBRA_BLOCKDIAGMATRIX              (BlockDiagMatrix)          " "${TEST_ALGEBRA_BLOCKDIAGMATRIX}")
printinfo_onoff("  TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX       (BlockGraphSparseMatrix)   " "${TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX}")
#printinfo_onoff("  TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX       (BlockLaplaceFreeMatrix)   " "${TEST_ALGEBRA_BLOCKLAPLACEFREEMATRIX}")
#printinfo_onoff("  TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX     (BlockLaplaceSparseMatrix) " "${TEST_ALGEBRA_BLOCKLAPLACESPARSEMATRIX}")
#printinfo_onoff("  TEST_ALGEBRA_DECOMPOSITION                (Decomposition)            " "${TEST_ALGEBRA_DECOMPOSITION}")
#printinfo_onoff("  TEST_ALGEBRA_FILECRSMATRIX                (FileCRSMatrix)            " "${TEST_ALGEBRA_FILECRSMATRIX}")
#printinfo_onoff("  TEST_ALGEBRA_PETSCVECTOR                  (PetscVector)              " "${TEST_ALGEBRA_PETSCVECTOR}")
#printinfo_onoff("  TEST_ALGEBRA_LOCALDENSEMATRIX             (LocalDenseMatrix)         " "${TEST_ALGEBRA_LOCALDENSEMATRIX}")
#printinfo_onoff("  TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL      (SimplexFeasibleSet_Local) " "${TEST_ALGEBRA_SIMPLEXFEASIBLESETLOCAL}")
printinfo_onoff(" TEST_DATA                                 (...)                      " "${TEST_DATA}")
#printinfo_onoff("  TEST_DATA_DIAG                            (DiagData)                 " "${TEST_DATA_DIAG}")
#printinfo_onoff("  TEST_DATA_EDF                             (EdfData)                  " "${TEST_DATA_EDF}")
#printinfo_onoff("  TEST_DATA_IMAGE                           (ImageData)                " "${TEST_DATA_IMAGE}")
#printinfo_onoff("  TEST_DATA_KMEANS                          (KmeansData)               " "${TEST_DATA_KMEANS}")
#printinfo_onoff("  TEST_DATA_QP                              (QPData)                   " "${TEST_DATA_QP}")
#printinfo_onoff("  TEST_DATA_SIGNAL1D                        (Signal1DData)             " "${TEST_DATA_SIGNAL1D}")
#printinfo_onoff("  TEST_DATA_SIMPLE                          (SimpleData)               " "${TEST_DATA_SIMPLE}")
#printinfo_onoff("  TEST_DATA_TS                              (TSData)                   " "${TEST_DATA_TS}")
printinfo_onoff(" TEST_MODEL                                (...)                      " "${TEST_MODEL}")
#printinfo_onoff("  TEST_MODEL_GRAPHH1FEM                     (GraphH1FEMModel)          " "${TEST_MODEL_GRAPHH1FEM}")
#printinfo_onoff("  TEST_MODEL_KMEANSH1FEM                    (KmeansH1FEMModel)         " "${TEST_MODEL_KMEANSH1FEM}")
printinfo_onoff(" TEST_SOLVER                               (...)                      " "${TEST_SOLVER}")
#printinfo_onoff("  TEST_SOLVER_CGQP                          (CGQPSolver)               " "${TEST_SOLVER_CGQP}")
#printinfo_onoff("  TEST_SOLVER_DIAG                          (DiagSolver)               " "${TEST_SOLVER_DIAG}")
#printinfo_onoff("  TEST_SOLVER_MULTICG                       (MultiCGSolver)            " "${TEST_SOLVER_MULTICG}")
#printinfo_onoff("  TEST_SOLVER_SIMPLE                        (SimpleSolver)             " "${TEST_SOLVER_SIMPLE}")
#printinfo_onoff("  TEST_SOLVER_SPGQP                         (SPGQPSolver)              " "${TEST_SOLVER_SPGQP}")
printinfo_onoff(" TEST_DLIB                                 (...)                      " "${TEST_DLIB}")
printinfo_onoff("  TEST_DLIB_ANNA                            (benchmark from Anna)       " "${TEST_DLIB_ANNA}")
printinfo_onoff("  TEST_DLIB_INTEGRAL                        (numerical integration)     " "${TEST_DLIB_INTEGRAL}")
printinfo_onoff("  TEST_DLIB_GUI                             (fun with X11)              " "${TEST_DLIB_GUI}")
 

# ----- COMMON -----

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

# ----- ALGEBRA -----

if(${TEST_ALGEBRA_DOT})
	# dot product
	if(${USE_CUDA})
		pascadd_executable("test_classes/algebra/test_dot.cu" "test_dot")
	else()
		pascadd_executable("test_classes/algebra/test_dot.cpp" "test_dot")
	endif()
	
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

if(${TEST_ALGEBRA_BLOCKGRAPHSPARSEMATRIX})
	# BlockGraphSparseMatrix
	if(${USE_CUDA})
		pascadd_executable("test_classes/algebra/test_blockgraphsparsematrix.cu" "test_blockgraphsparsematrix")
	else()
		pascadd_executable("test_classes/algebra/test_blockgraphsparsematrix.cpp" "test_blockgraphsparsematrix")
	endif()

	# copy data with sample graphs
	file(COPY "test_classes/data/test_algebra_blockgraphsparse/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

# ----- DATA -----

# ----- MODEL -----

# ----- SOLVER -----

# ----- DLIB ------
if(${TEST_DLIB_ANNA})
	# benchmark from anna - first experiences with dlib
	if(${USE_CUDA})
		pascadd_executable("test_classes/dlib/test_anna.cu" "test_anna")
	else()
		pascadd_executable("test_classes/dlib/test_anna.cpp" "test_anna")
	endif()

endif()

if(${TEST_DLIB_INTEGRAL})
	# test numerical integration
	if(${USE_CUDA})
		pascadd_executable("test_classes/dlib/test_integral.cu" "test_integral")
	else()
		pascadd_executable("test_classes/dlib/test_integral.cpp" "test_integral")
	endif()

endif()

if(${TEST_DLIB_GUI})
	# test graphical user interface
	if(${USE_CUDA})
		pascadd_executable("test_classes/dlib/test_gui.cu" "test_gui")
	else()
		pascadd_executable("test_classes/dlib/test_gui.cpp" "test_gui")
	endif()

endif()
