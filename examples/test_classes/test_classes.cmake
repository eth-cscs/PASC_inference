include_directories("${CMAKE_SOURCE_DIR}/test_classes/")

print("Classes tests")

#include seqarrayvector
include("test_classes/seqarrayvector/test_classes.cmake")

#include PETSC testing
if(${USE_PETSC})
	include("test_classes/petscvector/test_classes.cmake")
endif()
