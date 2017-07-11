include_directories("${CMAKE_SOURCE_DIR}/test_classes/petscvector/")

# ----- COMMON -----
option(TEST_PETSCVECTOR_COMMON 						"TEST_PETSCVECTOR_COMMON" OFF)
option(TEST_PETSCVECTOR_COMMON_CONSOLEARG 			  "TEST_PETSCVECTOR_COMMON_CONSOLEARG" OFF)
option(TEST_PETSCVECTOR_COMMON_CONSOLEOUTPUT		  "TEST_PETSCVECTOR_COMMON_CONSOLEOUTPUT" OFF)
option(TEST_PETSCVECTOR_COMMON_GLOBALMANAGER		  "TEST_PETSCVECTOR_COMMON_GLOBALMANAGER" OFF)
option(TEST_PETSCVECTOR_COMMON_LOGGING				  "TEST_PETSCVECTOR_COMMON_LOGGING" OFF)
option(TEST_PETSCVECTOR_COMMON_MEMORYCHECK			  "TEST_PETSCVECTOR_COMMON_MEMORYCHECK" OFF)
option(TEST_PETSCVECTOR_COMMON_OFFSET				  "TEST_PETSCVECTOR_COMMON_OFFSET" OFF)
option(TEST_PETSCVECTOR_COMMON_SHORTINFO			  "TEST_PETSCVECTOR_COMMON_SHORTINFO" OFF)
option(TEST_PETSCVECTOR_COMMON_TIMER				  "TEST_PETSCVECTOR_COMMON_TIMER" OFF)
if(${TEST_PETSCVECTOR_COMMON})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_COMMON_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- ALGEBRA -----
option(TEST_PETSCVECTOR_ALGEBRA						"TEST_PETSCVECTOR_ALGEBRA" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_DOT					  "TEST_PETSCVECTOR_ALGEBRA_DOT" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BGMGRAPH			  "TEST_PETSCVECTOR_ALGEBRA_BGMGRAPH" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID1D		  "TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID1D" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID2D		  "TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID2D" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BLOCKDIAGMATRIX		  "TEST_PETSCVECTOR_ALGEBRA_BLOCKDIAGMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX "TEST_PETSCVECTOR_ALGEBRA_BLOCKSPARSEMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX "TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX "TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_DECOMPOSITION		  "TEST_PETSCVECTOR_ALGEBRA_DECOMPOSITION" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_FILECRSMATRIX		  "TEST_PETSCVECTOR_ALGEBRA_FILECRSMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_LOCALDENSEMATRIX	  "TEST_PETSCVECTOR_ALGEBRA_LOCALDENSEMATRIX" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL "TEST_PETSCVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL" OFF)
option(TEST_PETSCVECTOR_ALGEBRA_FEM_IMAGE             "TEST_PETSCVECTOR_ALGEBRA_FEM_IMAGE" OFF)

if(${TEST_PETSCVECTOR_ALGEBRA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_ALGEBRA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- DATA -----
option(TEST_PETSCVECTOR_DATA 						"TEST_PETSCVECTOR_DATA" OFF)
option(TEST_PETSCVECTOR_DATA_DIAG					  "TEST_PETSCVECTOR_DATA_DIAG" OFF)
option(TEST_PETSCVECTOR_DATA_EDF					  "TEST_PETSCVECTOR_DATA_EDF" OFF)
option(TEST_PETSCVECTOR_DATA_IMAGE					  "TEST_PETSCVECTOR_DATA_IMAGE" OFF)
option(TEST_PETSCVECTOR_DATA_KMEANS					  "TEST_PETSCVECTOR_DATA_KMEANS" OFF)
option(TEST_PETSCVECTOR_DATA_QP						  "TEST_PETSCVECTOR_DATA_QP" OFF)
option(TEST_PETSCVECTOR_DATA_SIGNAL1D				  "TEST_PETSCVECTOR_DATA_SIGNAL1D" OFF)
option(TEST_PETSCVECTOR_DATA_SIMPLE					  "TEST_PETSCVECTOR_DATA_SIMPLE" OFF)
option(TEST_PETSCVECTOR_DATA_TS						  "TEST_PETSCVECTOR_DATA_TS" OFF)
if(${TEST_PETSCVECTOR_DATA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_DATA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- MODEL -----
option(TEST_PETSCVECTOR_MODEL						"TEST_PETSCVECTOR_MODEL" OFF)
option(TEST_PETSCVECTOR_MODEL_GRAPHH1FEM			  "TEST_PETSCVECTOR_MODEL_GRAPHH1FEM" OFF)
option(TEST_PETSCVECTOR_MODEL_KMEANSH1FEM			  "TEST_PETSCVECTOR_MODEL_KMEANSH1FEM" OFF)
if(${TEST_MODEL})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_MODEL_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- SOLVER -----
option(TEST_PETSCVECTOR_SOLVER						"TEST_PETSCVECTOR_SOLVER" OFF)
option(TEST_PETSCVECTOR_SOLVER_CGQP					  "TEST_PETSCVECTOR_SOLVER_CGQP" OFF)
option(TEST_PETSCVECTOR_SOLVER_DIAG					  "TEST_PETSCVECTOR_SOLVER_DIAG" OFF)
option(TEST_PETSCVECTOR_SOLVER_MULTICG				  "TEST_PETSCVECTOR_SOLVER_MULTICG" OFF)
option(TEST_PETSCVECTOR_SOLVER_SIMPLE				  "TEST_PETSCVECTOR_SOLVER_SIMPLE" OFF)
option(TEST_PETSCVECTOR_SOLVER_SPGQP				  "TEST_PETSCVECTOR_SOLVER_SPGQP" OFF)
if(${TEST_PETSCVECTOR_SOLVER})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_SOLVER_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- DLIB -----
option(TEST_PETSCVECTOR_DLIB						"TEST_PETSCVECTOR_DLIB" OFF)
option(TEST_PETSCVECTOR_DLIB_ANNA					  "TEST_PETSCVECTOR_DLIB_ANNA" OFF)
option(TEST_PETSCVECTOR_DLIB_INTEGRAL				  "TEST_PETSCVECTOR_DLIB_INTEGRAL" OFF)
option(TEST_PETSCVECTOR_DLIB_GUI					  "TEST_PETSCVECTOR_DLIB_GUI" OFF)
if(${TEST_PETSCVECTOR_DLIB})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_DLIB_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# ----- CUBA -----
option(TEST_PETSCVECTOR_CUBA						"TEST_PETSCVECTOR_CUBA" OFF)
option(TEST_PETSCVECTOR_CUBA_DEMO					  "TEST_PETSCVECTOR_CUBA_DEMO" OFF)
option(TEST_PETSCVECTOR_CUBA_ANNA					  "TEST_PETSCVECTOR_CUBA_ANNA" OFF)
if(${TEST_PETSCVECTOR_CUBA})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_CUBA_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# ----- general switching option
option(TEST_PETSCVECTOR "TEST_PETSCVECTOR" OFF)
if(${TEST_PETSCVECTOR})
	# define shortcut to compile all tests of this group
	getListOfVarsStartingWith("TEST_PETSCVECTOR_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()

# print info
printinfo_onoff(" TEST_PETSCVECTOR                                      (...)                          " "${TEST_PETSCVECTOR}")
printinfo_onoff("   TEST_PETSCVECTOR_COMMON                               (...)                        " "${TEST_PETSCVECTOR_COMMON}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_CONSOLEARG                    (ConsoleArg)               " "${TEST_PETSCVECTOR_COMMON_CONSOLEARG}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_CONSOLEOUTPUT                 (ConsoleOutput)            " "${TEST_PETSCVECTOR_COMMON_CONSOLEOUTPUT}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_GLOBALMANAGER                 (GlobalManager)            " "${TEST_PETSCVECTOR_COMMON_GLOBALMANAGER}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_LOGGING                       (Logging)                  " "${TEST_PETSCVECTOR_COMMON_LOGGING}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_MEMORYCHECK                   (MemoryCheck)              " "${TEST_PETSCVECTOR_COMMON_MEMORYCHECK}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_OFFSET                        (Offset)                   " "${TEST_PETSCVECTOR_COMMON_OFFSET}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_SHORTINFO                     (Shortinfo)                " "${TEST_PETSCVECTOR_COMMON_SHORTINFO}")
#printinfo_onoff("     TEST_PETSCVECTOR_COMMON_TIMER                         (Timer,StackTimer)         " "${TEST_PETSCVECTOR_COMMON_TIMER}")
printinfo_onoff("   TEST_PETSCVECTOR_ALGEBRA                              (...)                        " "${TEST_PETSCVECTOR_ALGEBRA}")
printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_DOT                          (dot product)              " "${TEST_PETSCVECTOR_ALGEBRA_DOT}")
printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_FEM_IMAGE                    (fem on images)            " "${TEST_PETSCVECTOR_ALGEBRA_FEM_IMAGE}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BGMGRAPH                     (BGMGraph)                 " "${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPH}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID1D               (BGMGraphGrid1D)           " "${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID1D}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID2D               (BGMGraphGrid2D)           " "${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID2D}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BLOCKDIAGMATRIX              (BlockDiagMatrix)          " "${TEST_PETSCVECTOR_ALGEBRA_BLOCKDIAGMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX       (BlockGraphSparseMatrix)   " "${TEST_PETSCVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX       (BlockLaplaceFreeMatrix)   " "${TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACEFREEMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX     (BlockLaplaceSparseMatrix) " "${TEST_PETSCVECTOR_ALGEBRA_BLOCKLAPLACESPARSEMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_DECOMPOSITION                (Decomposition)            " "${TEST_PETSCVECTOR_ALGEBRA_DECOMPOSITION}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_FILECRSMATRIX                (FileCRSMatrix)            " "${TEST_PETSCVECTOR_ALGEBRA_FILECRSMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_PETSCVECTOR                  (PetscVector)              " "${TEST_PETSCVECTOR_ALGEBRA_PETSCVECTOR}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_LOCALDENSEMATRIX             (LocalDenseMatrix)         " "${TEST_PETSCVECTOR_ALGEBRA_LOCALDENSEMATRIX}")
#printinfo_onoff("     TEST_PETSCVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL      (SimplexFeasibleSet_Local) " "${TEST_PETSCVECTOR_ALGEBRA_SIMPLEXFEASIBLESETLOCAL}")
printinfo_onoff("   TEST_PETSCVECTOR_DATA                                 (...)                        " "${TEST_PETSCVECTOR_DATA}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_DIAG                            (DiagData)                 " "${TEST_PETSCVECTOR_DATA_DIAG}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_EDF                             (EdfData)                  " "${TEST_PETSCVECTOR_DATA_EDF}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_IMAGE                           (ImageData)                " "${TEST_PETSCVECTOR_DATA_IMAGE}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_KMEANS                          (KmeansData)               " "${TEST_PETSCVECTOR_DATA_KMEANS}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_QP                              (QPData)                   " "${TEST_PETSCVECTOR_DATA_QP}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_SIGNAL1D                        (Signal1DData)             " "${TEST_PETSCVECTOR_DATA_SIGNAL1D}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_SIMPLE                          (SimpleData)               " "${TEST_PETSCVECTOR_DATA_SIMPLE}")
#printinfo_onoff("     TEST_PETSCVECTOR_DATA_TS                              (TSData)                   " "${TEST_PETSCVECTOR_DATA_TS}")
printinfo_onoff("   TEST_PETSCVECTOR_MODEL                                (...)                        " "${TEST_PETSCVECTOR_MODEL}")
#printinfo_onoff("     TEST_PETSCVECTOR_MODEL_GRAPHH1FEM                     (GraphH1FEMModel)          " "${TEST_PETSCVECTOR_MODEL_GRAPHH1FEM}")
#printinfo_onoff("     TEST_PETSCVECTOR_MODEL_KMEANSH1FEM                    (KmeansH1FEMModel)         " "${TEST_PETSCVECTOR_MODEL_KMEANSH1FEM}")
printinfo_onoff("   TEST_PETSCVECTOR_SOLVER                               (...)                        " "${TEST_PETSCVECTOR_SOLVER}")
#printinfo_onoff("     TEST_PETSCVECTOR_SOLVER_CGQP                          (CGQPSolver)               " "${TEST_PETSCVECTOR_SOLVER_CGQP}")
#printinfo_onoff("     TEST_PETSCVECTOR_SOLVER_DIAG                          (DiagSolver)               " "${TEST_PETSCVECTOR_SOLVER_DIAG}")
#printinfo_onoff("     TEST_PETSCVECTOR_SOLVER_MULTICG                       (MultiCGSolver)            " "${TEST_PETSCVECTOR_SOLVER_MULTICG}")
#printinfo_onoff("     TEST_PETSCVECTOR_SOLVER_SIMPLE                        (SimpleSolver)             " "${TEST_PETSCVECTOR_SOLVER_SIMPLE}")
#printinfo_onoff("     TEST_PETSCVECTOR_SOLVER_SPGQP                         (SPGQPSolver)              " "${TEST_PETSCVECTOR_SOLVER_SPGQP}")
printinfo_onoff("   TEST_PETSCVECTOR_DLIB                                 (...)                        " "${TEST_PETSCVECTOR_DLIB}")
#printinfo_onoff("     TEST_PETSCVECTOR_DLIB_ANNA                            (benchmark from Anna)      " "${TEST_PETSCVECTOR_DLIB_ANNA}")
#printinfo_onoff("     TEST_PETSCVECTOR_DLIB_INTEGRAL                        (numerical integration)    " "${TEST_PETSCVECTOR_DLIB_INTEGRAL}")
#printinfo_onoff("     TEST_PETSCVECTOR_DLIB_GUI                             (fun with X11)             " "${TEST_PETSCVECTOR_DLIB_GUI}")
printinfo_onoff("   TEST_PETSCVECTOR_CUBA                                 (...)                        " "${TEST_PETSCVECTOR_CUBA}")
printinfo_onoff("     TEST_PETSCVECTOR_CUBA_DEMO                            (demo from cuba library)    " "${TEST_PETSCVECTOR_CUBA_DEMO}")
printinfo_onoff("     TEST_PETSCVECTOR_CUBA_ANNA                            (benchmark from Anna)       " "${TEST_PETSCVECTOR_CUBA_ANNA}")
 

# ----- COMMON -----

if(${TEST_PETSCVECTOR_COMMON_CONSOLEARG})
	# ConsoleArgClass
	testadd_executable("test_classes/petscvector/common/test_consolearg.cpp" "test_petscvector_consolearg")
endif()

if(${TEST_PETSCVECTOR_COMMON_CONSOLEOUTPUT})
	# ConsoleOutput
	testadd_executable("test_classes/petscvector/common/test_consoleoutput.cpp" "test_petscvector_consoleoutput")
endif()

if(${TEST_PETSCVECTOR_COMMON_GLOBALMANAGER})
	# GlobalManager
	testadd_executable("test_classes/petscvector/common/test_globalmanager.cpp" "test_petscvector_globalmanager")
endif()

if(${TEST_PETSCVECTOR_COMMON_LOGGING})
	# Logging
	testadd_executable("test_classes/petscvector/common/test_logging.cpp" "test_petscvector_logging")
endif()

if(${TEST_PETSCVECTOR_COMMON_MEMORYCHECK})
	# MemoryCheck
	testadd_executable("test_classes/petscvector/common/test_memorycheck.cpp" "test_petscvector_memorycheck")
endif()

if(${TEST_PETSCVECTOR_COMMON_OFFSET})
	# Offset
	testadd_executable("test_classes/petscvector/common/test_offset.cpp" "test_petscvector_offset")
endif()

if(${TEST_PETSCVECTOR_COMMON_SHORTINFO})
	# Shortinfo
	testadd_executable("test_classes/petscvector/common/test_shortinfo.cpp" "test_petscvector_shortinfo")
endif()

if(${TEST_PETSCVECTOR_COMMON_TIMER})
	# Timer and StackTimer
	testadd_executable("test_classes/petscvector/common/test_timer.cpp" "test_petscvector_timer")
endif()

# ----- ALGEBRA -----

if(${TEST_PETSCVECTOR_ALGEBRA_DOT})
	# dot product
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/algebra/test_dot.cu" "test_petscvector_dot")
	else()
		testadd_executable("test_classes/petscvector/algebra/test_dot.cpp" "test_petscvector_dot")
	endif()
	
endif()

if(${TEST_PETSCVECTOR_ALGEBRA_FEM_IMAGE})
	# dot product
	testadd_executable("test_classes/petscvector/algebra/test_fem_image.cpp" "test_petscvector_fem_image")

	# copy scripts
	make_directory("scripts/test_classes/")
	file(COPY "scripts/" DESTINATION "scripts/test_classes/"	FILES_MATCHING PATTERN "*")
	file(COPY "test_classes/petscvector/scripts/" DESTINATION "scripts/test_classes/" FILES_MATCHING PATTERN "*")
	
	# copy data from test_image
	make_directory("data/test_image/")
	file(COPY "test_image/data/" DESTINATION "data/test_image/" FILES_MATCHING PATTERN "*")
endif()

if(${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPH})
	# BGMGraph
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraph.cu" "test_petscvector_bgmgraph")
	else()
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraph.cpp" "test_petscvector_bgmgraph")
	endif()
	
	# copy data with sample graphs
	file(COPY "test_classes/petscvector/data/test_algebra_bgmgraph/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

if(${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID1D})
	# BGMGraphGrid1D
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraphgrid1D.cu" "test_petscvector_bgmgraphgrid1D")
	else()
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraphgrid1D.cpp" "test_petscvector_bgmgraphgrid1D")
	endif()
endif()

if(${TEST_PETSCVECTOR_ALGEBRA_BGMGRAPHGRID2D})
	# BGMGraphGrid2D
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraphgrid2D.cu" "test_petscvector_bgmgraphgrid2D")
	else()
		testadd_executable("test_classes/petscvector/algebra/test_bgmgraphgrid2D.cpp" "test_petscvector_bgmgraphgrid2D")
	endif()
endif()

if(${TEST_PETSCVECTOR_ALGEBRA_BLOCKGRAPHSPARSEMATRIX})
	# BlockGraphSparseMatrix
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/algebra/test_blockgraphsparsematrix.cu" "test_petscvector_blockgraphsparsematrix")
	else()
		testadd_executable("test_classes/petscvector/algebra/test_blockgraphsparsematrix.cpp" "test_petscvector_blockgraphsparsematrix")
	endif()

	# copy data with sample graphs
	file(COPY "test_classes/petscvector/data/test_algebra_blockgraphsparse/" DESTINATION "data" FILES_MATCHING PATTERN "*")
	
endif()

# ----- DATA -----

# ----- MODEL -----

# ----- SOLVER -----

# ----- DLIB ------
if(${TEST_PETSCVECTOR_DLIB_ANNA})
	# benchmark from anna - first experiences with dlib
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/dlib/test_anna.cu" "test_petscvector_anna")
	else()
		testadd_executable("test_classes/petscvector/dlib/test_anna.cpp" "test_petscvector_anna")
	endif()

endif()

if(${TEST_PETSCVECTOR_DLIB_INTEGRAL})
	# test numerical integration
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/dlib/test_integral.cu" "test_petscvector_integral")
	else()
		testadd_executable("test_classes/petscvector/dlib/test_integral.cpp" "test_petscvector_integral")
	endif()

endif()

if(${TEST_PETSCVECTOR_DLIB_GUI})
	# test graphical user interface
	if(${USE_CUDA})
		testadd_executable("test_classes/petscvector/dlib/test_gui.cu" "test_petscvector_gui")
	else()
		testadd_executable("test_classes/petscvector/dlib/test_gui.cpp" "test_petscvector_gui")
	endif()

endif()

# ----- CUBA ------
if(${TEST_PETSCVECTOR_CUBA_DEMO})
	testadd_executable("test_classes/petscvector/cuba/test_cuba.cpp" "test_petscvector_cuba_demo")
endif()

if(${TEST_PETSCVECTOR_CUBA_ANNA})
	include_directories("${CMAKE_SOURCE_DIR}/test_classes/petscvector/cuba/include")

	add_executable(test_petscvector_cuba_anna test_classes/petscvector/cuba/test_anna.cpp test_classes/petscvector/cuba/ExtraParameters.cpp test_classes/petscvector/cuba/Integrator.cpp)
	set_source_files_properties("test_classes/petscvector/cuba/test_anna.cpp" COMPILE_FLAGS "${FLAGS_DEF_D}")
	set_target_properties(test_petscvector_cuba_anna PROPERTIES	OUTPUT_NAME "test_petscvector_cuba_anna")
	target_link_libraries(test_petscvector_cuba_anna ${LIBRARIES_DEF})

endif()

