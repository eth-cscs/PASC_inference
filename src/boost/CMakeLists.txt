project(PASC_INFERENCE_LIBRARY_BOOST)
cmake_minimum_required(VERSION 2.8.11)

# compile cpp files with boost 
if(${USE_BOOST})
	file(GLOB_RECURSE SOURCES_BOOST "${PASCINFERENCE_SRC}/boost/*.cpp")
	set_source_files_properties(${SOURCES_BOOST} COMPILE_FLAGS "${FLAGS_DEF_D}")
	add_library(pascinference_boost SHARED ${SOURCES_BOOST})

endif()

