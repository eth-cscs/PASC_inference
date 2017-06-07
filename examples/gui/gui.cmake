include_directories("${CMAKE_SOURCE_DIR}/gui/")

# decide which gui util to compile

option(GUI "GUI" OFF)
option(GUI_SHOW_VEC "GUI_SHOW_VEC" OFF)
option(GUI_SHOW_GAMMA "GUI_SHOW_GAMMA" OFF)
option(GUI_SHOW_IMAGE "GUI_SHOW_IMAGE" OFF)
option(GUI_SHOW_EEG "GUI_SHOW_NEURO" OFF)
option(GUI_SHOW_SHORTINFO "GUI_SHOW_SHORTINFO" OFF)

if(${GUI})
	# define shortcut to compile all utils
	getListOfVarsStartingWith("GUI_" matchedVars)
	foreach (_var IN LISTS matchedVars)
		set(${_var} ON)
	endforeach()
endif()


# print info
print("Graphical User Interface")
printinfo_onoff(" GUI                                                   (...)                          " "${GUI}")
printinfo_onoff("   GUI_SHOW_VEC                                          (ConsoleArg)                 " "${GUI_SHOW_VEC}")
printinfo_onoff("   GUI_SHOW_GAMMA                                        (ConsoleArg)                 " "${GUI_SHOW_GAMMA}")
printinfo_onoff("   GUI_SHOW_IMAGE                                        (ConsoleArg)                 " "${GUI_SHOW_IMAGE}")
printinfo_onoff("   GUI_SHOW_NEURO                                        (ConsoleArg)                 " "${GUI_SHOW_NEURO}")
printinfo_onoff("   GUI_SHOW_SHORTINFO                                    (ConsoleArg)                 " "${GUI_SHOW_SHORTINFO}")

if(${GUI_SHOW_VEC})
	set(CMAKE_BUILD_TYPE Release)
	testadd_executable("gui/gui_show_vec.cpp" "gui_show_vec")
endif()

if(${GUI_SHOW_GAMMA})
	set(CMAKE_BUILD_TYPE Release)
	testadd_executable("gui/gui_show_gamma.cpp" "gui_show_gamma")
endif()

if(${GUI_SHOW_IMAGE})
	set(CMAKE_BUILD_TYPE Release)
	testadd_executable("gui/gui_show_image.cpp" "gui_show_image")
endif()

if(${GUI_SHOW_NEURO})
	set(CMAKE_BUILD_TYPE Release)
	testadd_executable("gui/gui_show_neuro.cpp" "gui_show_neuro")
endif()

if(${GUI_SHOW_SHORTINFO})
	set(CMAKE_BUILD_TYPE Release)
	testadd_executable("gui/gui_show_shortinfo.cpp" "gui_show_shortinfo")
endif()
