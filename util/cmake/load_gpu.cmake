
if(${USE_GPU})

	#we have to use CUDA to run on GPU
	if(NOT ${USE_CUDA})
		message(FATAL_ERROR "${Red}Sorry, you cannot use USE_GPU without CUDA! (use -DCUDA=ON)${ColourReset}")
	endif()

	set(FLAGS_DEF "-USE_GPU ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_GPU ${FLAGS_DEF_D}")

endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_GPU)
	printinfo_onoff("USE_GPU\t\t\t" "${USE_GPU}")
endmacro()

