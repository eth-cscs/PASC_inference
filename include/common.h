
#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

/* general default values */
#define DEFAULT_DEBUG_MODE 0

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export results to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */


/* include MINLIN */ 
#ifdef USE_MINLIN
	#include <minlin/minlin.h>
	#include <minlin/modules/threx/threx.h>

	namespace minlin {
		namespace threx {
			MINLIN_INIT
		}
	}
#endif

/* PetscVector */
#ifdef USE_PETSC
	#include "petsc.h"
	#include "petscvector.h"

#endif


/* cuda stuff */
#ifdef USE_GPU
	#include <stdio.h> /* printf in cuda */

	#include <cuda.h>
	#include <cuda_runtime_api.h>
	#include <device_launch_parameters.h>
	#include <device_functions.h>
#endif

/* we are using namespace pascinference */
namespace pascinference {

/* global variables */
int DEBUG_MODE;

/* general utils */
void Initialize(int, char**);
void Finalize();

void Message(std::string text);
void Message_info(std::string text);
void Message_info_main(std::string text);
void Message_info_values(std::string text, std::string values);
void Message_info_value(std::string text, int value);
void Message_info_value(std::string text, double value);
void Message_info_time(std::string text, double value);
//void Message_error(string text);


/* structure for time management (measure computing time) */
class StackTimer {
		std::stack<double> time_stack;
		double time;

		double getUnixTime(void);
	public:
		void start();
		double stop();
		int status();
	
};

class Timer {
		double time_sum;
		double time_start;
		double time_last;

		double getUnixTime(void);
		bool run_or_not;
	public:
		void restart();
		void start();
		void stop();
		double get_value_sum();
		double get_value_last();
		bool status();
};


#ifdef USE_GPU
	/* cuda error check */ 
	#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
	inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
	{
		if (code != cudaSuccess) 
		{
			fprintf(stderr,"\n\x1B[31mCUDA error:\x1B[0m %s %s \x1B[33m%d\x1B[0m\n\n", cudaGetErrorString(code), file, line);
			if (abort) exit(code);
		}
	}
#endif

}

#endif
