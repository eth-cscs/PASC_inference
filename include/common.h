/*! This class should include all common functions and includes */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

/* include common c++ header files */
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stack>
#include <limits>

/* include MINLIN */
#include <minlin/minlin.h>
#include <minlin/modules/threx/threx.h>
//#include <qpopt/smalbe.h>

#include <stdio.h> /* printf in cuda */

/* cuda stuff */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <device_functions.h>

/* include settings */
#include "settings.h"

using namespace minlin::threx;
MINLIN_INIT


/* for global time management */
class Timer;
extern Timer timer;

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
class Timer {
		std::stack<double> time_stack;
		double time;

		double getUnixTime(void);
	public:
		void start();
		double stop();
		int status();
	
};

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
