include_directories("${CMAKE_SOURCE_DIR}/test_classes/seqarrayvector/")

# ----- COMMON -----
#option(TEST_SEQARRAYVECTOR_COMMON 						"TEST_SEQARRAYVECTOR_COMMON" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_CONSOLEARG 			  "TEST_SEQARRAYVECTOR_COMMON_CONSOLEARG" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_CONSOLEOUTPUT		  "TEST_SEQARRAYVECTOR_COMMON_CONSOLEOUTPUT" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_GLOBALMANAGER		  "TEST_SEQARRAYVECTOR_COMMON_GLOBALMANAGER" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_LOGGING				  "TEST_SEQARRAYVECTOR_COMMON_LOGGING" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_MEMORYCHECK			  "TEST_SEQARRAYVECTOR_COMMON_MEMORYCHECK" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_OFFSET				  "TEST_SEQARRAYVECTOR_COMMON_OFFSET" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_SHORTINFO			  "TEST_SEQARRAYVECTOR_COMMON_SHORTINFO" OFF)
#option(TEST_SEQARRAYVECTOR_COMMON_TIMER				  "TEST_SEQARRAYVECTOR_COMMON_TIMER" OFF)
if(${TEST_SEQARRAYVECTOR_COMMON})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_COMMON_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- ALGEBRA -----
option(TEST_SEQARRAYVECTOR_ALGEBRA						"TEST_SEQARRAYVECTOR_ALGEBRA" OFF)
option(TEST_SEQARRAYVECTOR_ALGEBRA_PRECISION					  "TEST_SEQARRAYVECTOR_ALGEBRA_PRECISION" OFF)
option(TEST_SEQARRAYVECTOR_ALGEBRA_DOT					  "TEST_SEQARRAYVECTOR_ALGEBRA_DOT" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPH			  "TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPH" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID1D		  "TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID1D" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID2D		  "TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID2D" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKDIAGMATRIX		  "TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKDIAGMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX "TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKSPARSEMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX "TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX "TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_DECOMPOSITION		  "TEST_SEQARRAYVECTOR_ALGEBRA_DECOMPOSITION" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_FILECRSMATRIX		  "TEST_SEQARRAYVECTOR_ALGEBRA_FILECRSMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_LOCALDENSEMATRIX	  "TEST_SEQARRAYVECTOR_ALGEBRA_LOCALDENSEMATRIX" OFF)
#option(TEST_SEQARRAYVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL "TEST_SEQARRAYVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL" OFF)
if(${TEST_SEQARRAYVECTOR_ALGEBRA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_ALGEBRA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- DATA -----
option(TEST_SEQARRAYVECTOR_DATA 						"TEST_SEQARRAYVECTOR_DATA" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_DIAG					  "TEST_SEQARRAYVECTOR_DATA_DIAG" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_EDF					  "TEST_SEQARRAYVECTOR_DATA_EDF" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_IMAGE					  "TEST_SEQARRAYVECTOR_DATA_IMAGE" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_KMEANS					  "TEST_SEQARRAYVECTOR_DATA_KMEANS" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_QP						  "TEST_SEQARRAYVECTOR_DATA_QP" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_SIGNAL1D				  "TEST_SEQARRAYVECTOR_DATA_SIGNAL1D" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_SIMPLE					  "TEST_SEQARRAYVECTOR_DATA_SIMPLE" OFF)
#option(TEST_SEQARRAYVECTOR_DATA_TS						  "TEST_SEQARRAYVECTOR_DATA_TS" OFF)
if(${TEST_SEQARRAYVECTOR_DATA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_DATA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- MODEL -----
option(TEST_SEQARRAYVECTOR_MODEL						"TEST_SEQARRAYVECTOR_MODEL" OFF)
#option(TEST_SEQARRAYVECTOR_MODEL_GRAPHH1FEM			  "TEST_SEQARRAYVECTOR_MODEL_GRAPHH1FEM" OFF)
#option(TEST_SEQARRAYVECTOR_MODEL_KMEANSH1FEM			  "TEST_SEQARRAYVECTOR_MODEL_KMEANSH1FEM" OFF)
if(${TEST_MODEL})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_MODEL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- SOLVER -----
option(TEST_SEQARRAYVECTOR_SOLVER						"TEST_SEQARRAYVECTOR_SOLVER" OFF)
#option(TEST_SEQARRAYVECTOR_SOLVER_CGQP					  "TEST_SEQARRAYVECTOR_SOLVER_CGQP" OFF)
#option(TEST_SEQARRAYVECTOR_SOLVER_DIAG					  "TEST_SEQARRAYVECTOR_SOLVER_DIAG" OFF)
#option(TEST_SEQARRAYVECTOR_SOLVER_MULTICG				  "TEST_SEQARRAYVECTOR_SOLVER_MULTICG" OFF)
#option(TEST_SEQARRAYVECTOR_SOLVER_SIMPLE				  "TEST_SEQARRAYVECTOR_SOLVER_SIMPLE" OFF)
#option(TEST_SEQARRAYVECTOR_SOLVER_SPGQP				  "TEST_SEQARRAYVECTOR_SOLVER_SPGQP" OFF)
if(${TEST_SEQARRAYVECTOR_SOLVER})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_SOLVER_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- DLIB -----
option(TEST_SEQARRAYVECTOR_DLIB						"TEST_SEQARRAYVECTOR_DLIB" OFF)
#option(TEST_SEQARRAYVECTOR_DLIB_ANNA					  "TEST_SEQARRAYVECTOR_DLIB_ANNA" OFF)
#option(TEST_SEQARRAYVECTOR_DLIB_INTEGRAL				  "TEST_SEQARRAYVECTOR_DLIB_INTEGRAL" OFF)
#option(TEST_SEQARRAYVECTOR_DLIB_GUI					  "TEST_SEQARRAYVECTOR_DLIB_GUI" OFF)
if(${TEST_SEQARRAYVECTOR_DLIB})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_DLIB_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- general switching option
option(TEST_SEQARRAYVECTOR "TEST_SEQARRAYVECTOR" OFF)
if(${TEST_SEQARRAYVECTOR})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_SEQARRAYVECTOR_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# print info
printinfo_onoff(" TEST_SEQARRAYVECTOR                                   (...)                          " "${TEST_SEQARRAYVECTOR}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_COMMON                            (...)                        " "${TEST_SEQARRAYVECTOR_COMMON}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_CONSOLEARG                 (ConsoleArg)               " "${TEST_SEQARRAYVECTOR_COMMON_CONSOLEARG}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_CONSOLEOUTPUT              (ConsoleOutput)            " "${TEST_SEQARRAYVECTOR_COMMON_CONSOLEOUTPUT}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_GLOBALMANAGER              (GlobalManager)            " "${TEST_SEQARRAYVECTOR_COMMON_GLOBALMANAGER}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_LOGGING                    (Logging)                  " "${TEST_SEQARRAYVECTOR_COMMON_LOGGING}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_MEMORYCHECK                (MemoryCheck)              " "${TEST_SEQARRAYVECTOR_COMMON_MEMORYCHECK}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_OFFSET                     (Offset)                   " "${TEST_SEQARRAYVECTOR_COMMON_OFFSET}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_SHORTINFO                  (Shortinfo)                " "${TEST_SEQARRAYVECTOR_COMMON_SHORTINFO}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_COMMON_TIMER                      (Timer,StackTimer)         " "${TEST_SEQARRAYVECTOR_COMMON_TIMER}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_ALGEBRA                           (...)                        " "${TEST_SEQARRAYVECTOR_ALGEBRA}")
printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_PRECISION                 (test precision)           " "${TEST_SEQARRAYVECTOR_ALGEBRA_PRECISION}")
printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_DOT                       (dot product)              " "${TEST_SEQARRAYVECTOR_ALGEBRA_DOT}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPH                  (BGMGraph)                 " "${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPH}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID1D            (BGMGraphGrid1D)           " "${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID1D}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID2D            (BGMGraphGrid2D)           " "${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID2D}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKDIAGMATRIX           (BlockDiagMatrix)          " "${TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKDIAGMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX    (BlockGraphSparseMatrix)   " "${TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX    (BlockLaplaceFreeMatrix)   " "${TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX  (BlockLaplaceSparseMatrix) " "${TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_DECOMPOSITION             (Decomposition)            " "${TEST_SEQARRAYVECTOR_ALGEBRA_DECOMPOSITION}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_FILECRSMATRIX             (FileCRSMatrix)            " "${TEST_SEQARRAYVECTOR_ALGEBRA_FILECRSMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_SEQARRAYVECTOR            (PetscVector)              " "${TEST_SEQARRAYVECTOR_ALGEBRA_SEQARRAYVECTOR}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_LOCALDENSEMATRIX          (LocalDenseMatrix)         " "${TEST_SEQARRAYVECTOR_ALGEBRA_LOCALDENSEMATRIX}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL   (SimplexFeasibleSet_Local) " "${TEST_SEQARRAYVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_DATA                              (...)                        " "${TEST_SEQARRAYVECTOR_DATA}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_DIAG                         (DiagData)                 " "${TEST_SEQARRAYVECTOR_DATA_DIAG}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_EDF                          (EdfData)                  " "${TEST_SEQARRAYVECTOR_DATA_EDF}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_IMAGE                        (ImageData)                " "${TEST_SEQARRAYVECTOR_DATA_IMAGE}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_KMEANS                       (KmeansData)               " "${TEST_SEQARRAYVECTOR_DATA_KMEANS}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_QP                           (QPData)                   " "${TEST_SEQARRAYVECTOR_DATA_QP}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_SIGNAL1D                     (Signal1DData)             " "${TEST_SEQARRAYVECTOR_DATA_SIGNAL1D}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_SIMPLE                       (SimpleData)               " "${TEST_SEQARRAYVECTOR_DATA_SIMPLE}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DATA_TS                           (TSData)                   " "${TEST_SEQARRAYVECTOR_DATA_TS}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_MODEL                             (...)                        " "${TEST_SEQARRAYVECTOR_MODEL}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_MODEL_GRAPHH1FEM                  (GraphH1FEMModel)          " "${TEST_SEQARRAYVECTOR_MODEL_GRAPHH1FEM}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_MODEL_KMEANSH1FEM                 (KmeansH1FEMModel)         " "${TEST_SEQARRAYVECTOR_MODEL_KMEANSH1FEM}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_SOLVER                            (...)                        " "${TEST_SEQARRAYVECTOR_SOLVER}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_SOLVER_CGQP                       (CGQPSolver)               " "${TEST_SEQARRAYVECTOR_SOLVER_CGQP}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_SOLVER_DIAG                       (DiagSolver)               " "${TEST_SEQARRAYVECTOR_SOLVER_DIAG}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_SOLVER_MULTICG                    (MultiCGSolver)            " "${TEST_SEQARRAYVECTOR_SOLVER_MULTICG}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_SOLVER_SIMPLE                     (SimpleSolver)             " "${TEST_SEQARRAYVECTOR_SOLVER_SIMPLE}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_SOLVER_SPGQP                      (SPGQPSolver)              " "${TEST_SEQARRAYVECTOR_SOLVER_SPGQP}")
printinfo_onoff("   TEST_SEQARRAYVECTOR_DLIB                              (...)                        " "${TEST_SEQARRAYVECTOR_DLIB}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DLIB_ANNA                         (benchmark from Anna)      " "${TEST_SEQARRAYVECTOR_DLIB_ANNA}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DLIB_INTEGRAL                     (numerical integration)    " "${TEST_SEQARRAYVECTOR_DLIB_INTEGRAL}")
#printinfo_onoff("     TEST_SEQARRAYVECTOR_DLIB_GUI                          (fun with X11)             " "${TEST_SEQARRAYVECTOR_DLIB_GUI}")
 

# ----- COMMON -----

if(${TEST_SEQARRAYVECTOR_COMMON_CONSOLEARG})
	# ConsoleArgClass
	testadd_executable("test_classes/seqarrayvector/common/test_consolearg.cpp" "test_seqarrayvector_consolearg")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_CONSOLEOUTPUT})
	# ConsoleOutput
	testadd_executable("test_classes/seqarrayvector/common/test_consoleoutput.cpp" "test_seqarrayvector_consoleoutput")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_GLOBALMANAGER})
	# GlobalManager
	testadd_executable("test_classes/seqarrayvector/common/test_globalmanager.cpp" "test_seqarrayvector_globalmanager")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_LOGGING})
	# Logging
	testadd_executable("test_classes/seqarrayvector/common/test_logging.cpp" "test_seqarrayvector_logging")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_MEMORYCHECK})
	# MemoryCheck
	testadd_executable("test_classes/seqarrayvector/common/test_memorycheck.cpp" "test_seqarrayvector_memorycheck")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_OFFSET})
	# Offset
	testadd_executable("test_classes/seqarrayvector/common/test_offset.cpp" "test_seqarrayvector_offset")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_SHORTINFO})
	# Shortinfo
	testadd_executable("test_classes/seqarrayvector/common/test_shortinfo.cpp" "test_seqarrayvector_shortinfo")
endif()

if(${TEST_SEQARRAYVECTOR_COMMON_TIMER})
	# Timer and StackTimer
	testadd_executable("test_classes/seqarrayvector/common/test_timer.cpp" "test_seqarrayvector_timer")
endif()

# ----- ALGEBRA -----

if(${TEST_SEQARRAYVECTOR_ALGEBRA_PRECISION})
	# BGMGraph
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/algebra/test_precision.cu" "test_seqarrayvector_precision")
	else()
		testadd_executable("test_classes/seqarrayvector/algebra/test_precision.cpp" "test_seqarrayvector_precision")
	endif()
endif()

if(${TEST_SEQARRAYVECTOR_ALGEBRA_DOT})
	# dot product
	testadd_executable("test_classes/seqarrayvector/algebra/test_dot.cpp" "test_seqarrayvector_dot")
endif()

if(${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPH})
	# BGMGraph
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraph.cu" "test_seqarrayvector_bgmgraph")
	else()
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraph.cpp" "test_seqarrayvector_bgmgraph")
	endif()
	
	# copy data with sample graphs
#	file(COPY "test_classes/seqarrayvector/data/test_algebra_bgmgraph/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

if(${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID1D})
	# BGMGraphGrid1D
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraphgrid1D.cu" "test_seqarrayvector_bgmgraphgrid1D")
	else()
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraphgrid1D.cpp" "test_seqarrayvector_bgmgraphgrid1D")
	endif()
endif()

if(${TEST_SEQARRAYVECTOR_ALGEBRA_BGMGRAPHGRID2D})
	# BGMGraphGrid2D
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraphgrid2D.cu" "test_seqarrayvector_bgmgraphgrid2D")
	else()
		testadd_executable("test_classes/seqarrayvector/algebra/test_bgmgraphgrid2D.cpp" "test_seqarrayvector_bgmgraphgrid2D")
	endif()
endif()

if(${TEST_SEQARRAYVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX})
	# BlockGraphSparseMatrix
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/algebra/test_blockgraphsparsematrix.cu" "test_seqarrayvector_blockgraphsparsematrix")
	else()
		testadd_executable("test_classes/seqarrayvector/algebra/test_blockgraphsparsematrix.cpp" "test_seqarrayvector_blockgraphsparsematrix")
	endif()

	# copy data with sample graphs
	file(COPY "test_classes/seqarrayvector/data/test_algebra_blockgraphsparse/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

# ----- DATA -----

# ----- MODEL -----

# ----- SOLVER -----

# ----- DLIB ------
if(${TEST_SEQARRAYVECTOR_DLIB_ANNA})
	# benchmark from anna - first experiences with dlib
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/dlib/test_anna.cu" "test_seqarrayvector_anna")
	else()
		testadd_executable("test_classes/seqarrayvector/dlib/test_anna.cpp" "test_seqarrayvector_anna")
	endif()

endif()

if(${TEST_SEQARRAYVECTOR_DLIB_INTEGRAL})
	# test numerical integration
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/dlib/test_integral.cu" "test_seqarrayvector_integral")
	else()
		testadd_executable("test_classes/seqarrayvector/dlib/test_integral.cpp" "test_seqarrayvector_integral")
	endif()

endif()

if(${TEST_SEQARRAYVECTOR_DLIB_GUI})
	# test graphical user interface
	if(${USE_CUDA})
		testadd_executable("test_classes/seqarrayvector/dlib/test_gui.cu" "test_seqarrayvector_gui")
	else()
		testadd_executable("test_classes/seqarrayvector/dlib/test_gui.cpp" "test_seqarrayvector_gui")
	endif()

endif()
