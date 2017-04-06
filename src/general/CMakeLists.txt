# compile general cpp files
file(GLOB SOURCES_GENERAL "${PASCINFERENCE_SRC}/general/*.cpp") # locate all cpp files in source directory
set_source_files_properties(${SOURCES_GENERAL} COMPILE_FLAGS ${FLAGS_DEF_D})
add_library(pascinference SHARED ${SOURCES_GENERAL})
