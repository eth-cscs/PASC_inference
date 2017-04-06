project(PASC_INFERENCE_LIBRARY_CUDA)
cmake_minimum_required(VERSION 2.8.11)

# compile cuda cu files
if(${USE_CUDA})
	file(GLOB_RECURSE SOURCES_CUDA "${PASCINFERENCE_SRC}/cuda/*.cu") # locate all cu files in source directory
	set_source_files_properties(${SOURCES_CUDA} COMPILE_FLAGS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\"")
	cuda_add_library(${PROJECT_NAME} OBJECT ${SOURCES_CUDA})
endif()

