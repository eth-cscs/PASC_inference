project(PASC_INFERENCE_LIBRARY_GENERAL)
cmake_minimum_required(VERSION 2.8.11)

# compile general cpp files
file(GLOB_RECURSE SOURCES_GENERAL "${PASCINFERENCE_SRC}/general/*.cpp") # locate all cpp files in source directory
set_source_files_properties(${SOURCES_GENERAL} COMPILE_FLAGS "${FLAGS_DEF_D}")
add_library(${PROJECT_NAME} OBJECT ${SOURCES_GENERAL})
