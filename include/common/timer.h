/** @file timer.h
 *  @brief for measuring time
 *
 *  This file includes classes and functions for manipulation with timers (StackTimer and Timer).
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_TIMER_H
#define	PASC_COMMON_TIMER_H

namespace pascinference {
namespace common {

/** \class StackTimer
 *  \brief stack-based time management
 *
 *  Class includes the stack with time markers. If user call start, then
 *  new record with the actual time is pushed into the stack. After calling
 *  stop, the difference between actual time and the top of stack is returned.
 *  Moreover, the top marker is removed.
 * 
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

		bool run_or_not; /**< is the timer running? */
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

}
} /* end of namespace */


/* ------------ IMPLEMENTATION ------------- */

namespace pascinference {
namespace common {
	
/* ------------ STACK TIMER ------------ */
double StackTimer::getUnixTime(void){
//	struct timespec tv;
//	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
//	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
	return MPI_Wtime();
}

void StackTimer::start(){
	double t_start = this->getUnixTime();
	this->time_stack.push(t_start);
}
	
double StackTimer::stop(){
	double t_end = this->getUnixTime();
	double t_start;
	if(this->time_stack.empty()){
		t_start = 0.0;
	} else {
		t_start = this->time_stack.top();
		this->time_stack.pop();
	}
	double out_time = double(t_end - t_start);
	
	return out_time;
}

int StackTimer::status(){
	return this->time_stack.size();
}

/* ------------ SIMPLE TIMER ------------ */
double Timer::getUnixTime(void){
//	struct timespec tv;
//	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
//	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
	return MPI_Wtime();
}

void Timer::restart(){
	this->time_sum = 0.0;
	this->time_last = 0.0;
	this->run_or_not = false;
	this->time_start = std::numeric_limits<double>::max();
}

void Timer::start(){
	this->time_start = this->getUnixTime();
	this->run_or_not = true;
}

void Timer::stop(){
	this->time_last = this->getUnixTime() - this->time_start;
	this->time_sum += this->time_last;
	this->run_or_not = false;
	this->time_start = std::numeric_limits<double>::max();
}

double Timer::get_value_sum() const {
	return this->time_sum;
}

double Timer::get_value_last() const {
	return this->time_last;
}

bool Timer::status() const {
	return this->run_or_not;
}


}
} /* end of namespace */


#endif
