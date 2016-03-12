/** @file common.h
 *  @brief commonly used stuff
 *
 *  This file includes commonly used classes and functions.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

/* general default values */
#define DEFAULT_DEBUG_MODE 0

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export solution to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */

#include "sys/types.h"
#include "sys/sysinfo.h"

/* we are using namespace pascinference */
namespace pascinference {

/* global variables */
int DEBUG_MODE; /**< the debug mode of the library */


/** @brief initialize the library
 * 
 *  Initialize random number generator. 
 *  If we are using Petsc, then Petsc is initialized.
 *
 * @todo process input arguments (use boost library?)
 */
void Initialize(int, char**); 

/** @brief finalize the library
 * 
 *  If we are using Petsc, then Petsc is finalized.
 *
 */
void Finalize();

void Message(std::string text);
void Message_info(std::string text);
void Message_info_main(std::string text);
void Message_info_values(std::string text, std::string values);
void Message_info_value(std::string text, int value);
void Message_info_value(std::string text, double value);
void Message_info_time(std::string text, double value);


/* structure for time management (measure computing time) */
/** \class StackTimer
 *  \brief stack-based time management
 *
 *  Class includes the stack with time markers. If user call start, then
 *  new record with the actual time is pushed into the stack. After calling
 *  stop, the difference between actual time and the top of stack is returned.
 *  Moreover, the top marker is removed.
*/
class StackTimer {
		std::stack<double> time_stack; /**< the stack with times */
		double time; 

		/** @brief get actual time
		* 
		*  Return the actual time in Unix-time format.
		*
		*/
		double getUnixTime(void); 

	public:

		/** @brief start to measure time
		* 
		*  Add actual time to the top of stack.
		*
		*/
		void start();

		/** @brief stop to measure time
		* 
		*  Return the difference between actual time and the top of the stack and
		*  remove the top node of the stack.
		*
		*/
		double stop();
		
		/** @brief get the size of the stack
		*/
		int status();
	
};

/** \class Timer
 *  \brief additive time management
 *
 *  Class measures the running time of the process. 
 *  It includes two time values - the total time and the time from the last call.
 *  User should use the timer as a sequences of start() and stop(). The last call 
 *  is the time between the last pair, total time is the time between all pairs.
 * 
*/
class Timer {
		double time_sum;
		double time_start;
		double time_last;

		/** @brief get actual time
		* 
		*  Return the actual time in Unix-time format.
		*
		*/
		double getUnixTime(void);

		bool run_or_not;
	public:

		/** @brief restart the timer
		* 
		*  Stop measuring time, set all time values equal to zero.
		*
		*/
		void restart();
		
		/** @brief start the timer
		* 
		*  Store actual time and start to measure time from this time.
		*
		*/
		void start();

		/** @brief store actual time
		* 
		*  Stop to measure time and compute the elapsed time.
		*  Store the ellapsed time as last time.
		*  Add elapsed time to the sum timer.
		* 
		*/
		void stop();
		
		/** @brief sum of timer calls
		* 
		*  Return the sum of ellapsed time between all start() and stop() calls.
		* 
		*/
		double get_value_sum();

		/** @brief get last elapsed time
		* 
		*  Return the last ellapsed time between start() and stop(). 
		* 
		*/
		double get_value_last();

		/** @brief get the status
		* 
		*  Return the status of timer 
		*  - true if the timer is running, 
		* false if the timer is not running.
		*
		*/
		bool status();
};


class MemoryCheck {
	public:
		static long long get_virtual_all() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalVirtualMem = memInfo.totalram;
			totalVirtualMem += memInfo.totalswap;
			totalVirtualMem *= memInfo.mem_unit;
			
			return totalVirtualMem;
		}

		static long long get_virtual_used() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
		    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;

			virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
			virtualMemUsed *= memInfo.mem_unit;
			
			return virtualMemUsed;
		}

		static double get_virtual(){
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			
			return 100*((memInfo.totalram - memInfo.freeram) + (memInfo.totalram - memInfo.freeram))/(double)(memInfo.totalram + memInfo.totalram);
		}

		static long long get_physical() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalPhysMem = memInfo.totalram;
			totalPhysMem *= memInfo.mem_unit;
			
			return totalPhysMem;
		}

		static void test_temp(int idx) {
			//std::cout << "TEST_" << idx << ": \033[31m" << get_virtual() << "\033[0m%" << std::endl;
		}
		
		
    //Add other values in next statement to avoid int overflow on right hand side...


	
	
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
