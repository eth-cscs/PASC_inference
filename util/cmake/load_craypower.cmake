# measure power consumption on Cray machines (like Piz Daint)
# based on code provided by Ben Cumming: https://github.com/bcumming/cray-power.git

if(${USE_CRAYPOWER})
	set(FLAGS_DEF "-USE_CRAYPOWER ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_CRAYPOWER ${FLAGS_DEF_D}")
endif()

# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_CRAYPOWER)
	printinfo_onoff("USE_CRAYPOWER\t\t\t" "${USE_CRAYPOWER}")
endmacro()

