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
#include <iostream>
#include <iomanip> /* setw - formated cout output */

/* we are using namespace pascinference */
namespace pascinference {

/* global variables */
int DEBUG_MODE; /**< the debug mode of the library */

#ifdef USE_PETSCVECTOR
 extern bool PETSC_INITIALIZED;
#endif

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
		double get_value_sum() const;

		/** @brief get last elapsed time
		* 
		*  Return the last ellapsed time between start() and stop(). 
		* 
		*/
		double get_value_last() const;

		/** @brief get the status
		* 
		*  Return the status of timer 
		*  - true if the timer is running, 
		* false if the timer is not running.
		*
		*/
		bool status() const;
};

/** @class MemoryCheck
 *  @brief get memory state
 * 
 *  Several utils for memory management. Could be used to conrol the memory.
 * 
 */ 
class MemoryCheck {
	public:
		/** @brief get the size of virtual memory
		 * 
		 */ 
		static long long get_virtual_all() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalVirtualMem = memInfo.totalram;
			totalVirtualMem += memInfo.totalswap;
			totalVirtualMem *= memInfo.mem_unit;
			
			return totalVirtualMem;
		}

		/** @brief get the size of used virtual memory
		 * 
		 */ 
		static long long get_virtual_used() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
		    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;

			virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
			virtualMemUsed *= memInfo.mem_unit;
			
			return virtualMemUsed;
		}

		/** @brief get the percentage usege of virtual memory
		 * 
		 */ 
		static double get_virtual(){
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			
			return 100*((memInfo.totalram - memInfo.freeram) + (memInfo.totalram - memInfo.freeram))/(double)(memInfo.totalram + memInfo.totalram);
		}

		/** @brief get the size of physical memory
		 * 
		 */ 
		static long long get_physical() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalPhysMem = memInfo.totalram;
			totalPhysMem *= memInfo.mem_unit;
			
			return totalPhysMem;
		}
	
};


class GlobalManagerClass {
	private:
		int rank;
		int size;
		bool initialized; /**< the rank and size were already obtained */

		void init(){
			if(!initialized){
			#ifdef USE_PETSCVECTOR
				/* can be set after initialize of petsc */
				if(petscvector::PETSC_INITIALIZED){
					TRY(PetscBarrier(NULL));
						
					MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
					MPI_Comm_size(MPI_COMM_WORLD, &this->size);

					initialized = true;
						
				}	
			#else
				initialized = true; /* if it is not with petsc, then this is always master */
				this->rank = 0;
				this->size = 1; /* one processor */
			#endif
			}
		}

	public:
		GlobalManagerClass(){
			this->init();
		}

		int get_rank(){
			this->init();
			return this->rank;
		}
		
		int get_size(){
			this->init();
			return this->size;
		}
};

static GlobalManagerClass GlobalManager; /**< for manipulation with rank and size of MPI */


/** @class OffsetClass
 *  @brief output space
 * 
 *  For manipulating with the space in the begining of output line.
 * 
 */ 
class OffsetClass {
	private:
		int size; /**< number of spaces in offset */

	public:
		/** @brief constructor
		* 
		*  Set initial size of offset to zero equal to zero.
		*
		*/
		OffsetClass() {
			size = 0;
		}
		
		/** @brief increase the size of offset
		*
		*/
		void push(){
			size += 3;
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			size -= 3;
			if(size < 0) size = 0;
		}

		friend std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset);

};

std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset){
	int i;
	for(i=0;i<my_offset.size;i++){
		output << " ";
	}
	return output;
}

OffsetClass offset; /**< global instance of output offset */


#ifdef USE_PETSC

/** @class ConsoleOutput
 *  @brief print only on master or on all processors
 * 
 *  Overloaded std::ostream cout function. Print based on rank.
 * 
 */ 
class ConsoleOutput : public std::ostream {
	private:

		/* Write a stream buffer that prefixes each line with Plop */
		class ConsoleOutputBuf: public std::stringbuf{
			private:
				std::ostream&   output;
			public:
				int print_rank; /**< which rank to print, -1 if all **/
				int rank; /**< rank of this process */

				ConsoleOutputBuf(std::ostream& str):output(str){
				}

				~ConsoleOutputBuf(){
				}

				virtual int sync ( ){
					if(this->rank == print_rank || this->print_rank == -1){
						std::stringstream output_string;

						/* write here also a rank of processor */
						if(GlobalManager.get_size() > 1){
							output_string << "[" << GlobalManager.get_rank() << "] " << offset << str();
						} else {
							output_string << offset << str();
						}
						str("");

//						output_string << output.rdbuf();

						if(this->print_rank == 0 || GlobalManager.get_size() <= 1){
							/* master prints */
							TRY( PetscPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
						} else {
							/* all prints */
							TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
						}

						output_string.str("");
						output_string.clear();

					} else {
					    str("");
					}

//					output.flush();
					return 0;
				}
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param std output stream (for example std::cout)
		*/
		ConsoleOutput(std::ostream& str, int rank = -1) : std::ostream(&buffer), buffer(str) {
			buffer.print_rank = rank;
			buffer.rank = GlobalManager.get_rank();

		}

		/** @brief increase the size of offset
		*
		*/
		void push(){
			offset.push();
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			offset.pop();
		}

		void synchronize(){
			if(buffer.print_rank == -1){
				TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
			}
		}

		
};

#else 

/** @class ConsoleOutput
 *  @brief print only on master or on all processors
 * 
 *  Overloaded std::ostream cout function. Print based on rank.
 * 
 */ 
class ConsoleOutput : public std::ostream {
	private:

		/* Write a stream buffer that prefixes each line with Rank */
		class ConsoleOutputBuf: public std::stringbuf{
			private:
				std::ostream&   output;
			public:
				int rank; /**< rank of this process */
				bool rankset; /**< the rank was already obtained */
				int print_rank; /**< which rank to print, -1 if all **/

				ConsoleOutputBuf(std::ostream& str):output(str){
				}

				~ConsoleOutputBuf(){
				}

				/** @brief set rank of this processor
				*
				*/
				void set_rank(){
					if(!rankset){
						#ifdef USE_PETSCVECTOR
							/* can be set after initialize of petsc */
							if(petscvector::PETSC_INITIALIZED){
								TRY(PetscBarrier(NULL));
						
								MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
								rankset = true;
						
							}	
						#else
							rankset = true; /* if it is not with petsc, then this is always master */
						#endif
					}
				}
				
				virtual int sync ( ){
					set_rank();
					if(this->rank == print_rank || this->print_rank == -1){
						#ifdef USE_PETSCVECTOR
							/* write here also a rank of processor */
							output << "[" << this->rank << "] " << offset << str();

						#else
							output << offset << str();
						#endif
						str("");
						output.flush();
					} else {
					    str("");
					}

					#ifdef USE_PETSCVECTOR
						TRY(PetscBarrier(NULL)); /* barrier after each cout */
					#endif

					return 0;
				}
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param std output stream (for example std::cout)
		*/
		ConsoleOutput(std::ostream& str, int rank = -1) : std::ostream(&buffer), buffer(str) {
			buffer.rankset = false;
			buffer.set_rank();
			buffer.print_rank = rank;
		}

		/** @brief increase the size of offset
		*
		*/
		void push(){
			offset.push();
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			offset.pop();
		}

		void synchronize(){
		}
		
};

#endif

static ConsoleOutput coutMaster(std::cout,0); /**< instance of output console stream on master */
static ConsoleOutput coutAll(std::cout); /**< instance of output console stream */


template<class MyType>
void print_array(std::ostream &output, int my_size, MyType *my_array){
	output << "[";
	for(int i=0;i<my_size;i++){
		output << my_array[i];
		if(i<my_size-1) output << ",";
	}	
	output << "]";
}

template<class MyType>
MyType sum_array(int my_size, const MyType *my_array){
	MyType sum = 0;
	for(int i=0;i<my_size;i++){
		sum += my_array[i];
	}	
	return sum;
}

template<class MyType>
MyType max_array(int my_size, const MyType *my_array){
	MyType return_value = my_array[0];
	for(int i=1;i<my_size;i++){
		if(my_array[i] > return_value){
			return_value = my_array[i];
		}
	}	
	return return_value;
}

template<class MyType>
MyType sum_subarray(int start, int end, const MyType *my_array){
	MyType sum = 0;
	for(int i=start;i<=end;i++){
		sum += my_array[i];
	}	
	return sum;
}

template<class MyType>
MyType dot_arrays(int size, const MyType *my_array1, const MyType *my_array2){
	MyType sum = 0;
	for(int i=0;i<=size;i++){
		sum += my_array1[i]*my_array2[i];
	}	
	return sum;
}

template<class MyType>
MyType max(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 > a1){
		return_value = a2;
	}
	return return_value;
}

template<class MyType>
MyType min(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 < a1){
		return_value = a2;
	}
	return return_value;
}


/* my_array[i] = my_array1[i]*my_array2[i] */
template<class MyType>
void mult_pw_array(int my_size, MyType *my_array, const MyType *my_array1, const MyType *my_array2){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_array1[i]*my_array2[i];
	}	
}

#ifdef USE_CUDA
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

#ifdef USE_PETSC
void MyVecDot(Vec x, Vec y, double *out){
	double *x_arr;
	double *y_arr;
	
	int size;
	TRY( VecGetSize(x, &size) );
	
	TRY( VecGetArray(x,&x_arr) );
	TRY( VecGetArray(y,&y_arr) );
	
	*out = dot_arrays(size, x_arr, y_arr);

	TRY( VecRestoreArray(x,&x_arr) );
	TRY( VecRestoreArray(y,&y_arr) );
}

/* w = x.*y */
void MyVecPointwiseMult(Vec w, Vec x, Vec y){
	double *w_arr;
	double *x_arr;
	double *y_arr;
	
	int size;
	TRY( VecGetSize(w, &size) );
	
	TRY( VecGetArray(w,&w_arr) );
	TRY( VecGetArray(x,&x_arr) );
	TRY( VecGetArray(y,&y_arr) );
	
	coutAll << "w: ";
	print_array(coutAll, size, w_arr);
	coutAll << std::endl;

	coutAll << "x: ";
	print_array(coutAll, size, x_arr);
	coutAll << std::endl;

	coutAll << "y: ";
	print_array(coutAll, size, x_arr);

	
	for(int i=0;i<=size;i++){
		w_arr[i] = x_arr[i]*y_arr[i];
	}	

	TRY( VecRestoreArray(w,&w_arr) );
	TRY( VecRestoreArray(x,&x_arr) );
	TRY( VecRestoreArray(y,&y_arr) );
}


#endif

}

#endif
