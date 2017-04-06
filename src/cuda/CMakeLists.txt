# compile cuda cu files
if($USE_CUDA)
	file(GLOB SOURCES_GENERAL "${PASCINFERENCE_SRC}/cuda/*.cu") # locate all cu files in source directory
	set_source_files_properties(${SOURCES_CUDA} COMPILE_FLAGS "${FLAGS_DEF_D} -arch=sm_60 --compiler-options \"${CUDA_CXX_FLAGS}\"")
	cuda_add_library(pascinference SHARED ${SOURCES_CUDA})
endif()

