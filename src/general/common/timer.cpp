#include "general/common/timer.h"

namespace pascinference {
namespace common {
	
/* ------------ STACK TIMER ------------ */
double StackTimer::getUnixTime(void){
//	struct timespec tv;
//	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
//	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
//	MPI_Barrier(MPI_COMM_WORLD);
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
//	MPI_Barrier(MPI_COMM_WORLD);
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

