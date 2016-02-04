#include "common.h"

/* global variables */
int DEBUG_MODE = DEFAULT_DEBUG_MODE; /* default debug mode */

/*!
 * initialize the application
 */ 
void Initialize(int argc, char *argv[]){

	/* initialize random seed: */
	if(RANDOM_BY_TIME){
		srand(time(NULL));
	} else {
		srand(0);
	}

  	/* init Petsc */
  	#ifdef USE_PETSC
		PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	#endif
}

/*!
 * final call of the application
 */ 
void Finalize(){
  	/* finalize Petsc */
  	#ifdef USE_PETSC
		PetscFinalize();
	#endif

}

/*!
 * print green message
 */
void Message(std::string text){
	std::cout << "\033[32m" << text << "\033[0m" << std::endl;
}  

/*!
 * print info message
 */
void Message_info_main(std::string text){
	std::cout << "\033[33m" << text << "\033[0m" << std::endl;
}  

void Message_info(std::string text){
	std::cout << "\033[36m" << text << "\033[0m" << std::endl;
}  

/*!
 * print info message with values and/or one value
 */
void Message_info_values(std::string text, std::string values){
	std::cout << "\033[36m" << text << "\033[0m" << values << std::endl;
}  
void Message_info_value(std::string text, int value){
	std::cout << "\033[36m" << text << "\033[0m" << value << std::endl;
}  
void Message_info_value(std::string text, double value){
	std::cout << "\033[36m" << text << "\033[0m" << value << std::endl;
}  
void Message_info_time(std::string text, double value){
	std::cout << "\033[36m" << text << "\033[0m" << value << "s" << std::endl;
}  





/* ------------ STACK TIMER ------------ */
double StackTimer::getUnixTime(void){
	struct timespec tv;
	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
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
	struct timespec tv;
	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
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

double Timer::get_value_sum(){
	return this->time_sum;
}

double Timer::get_value_last(){
	return this->time_last;
}

bool Timer::status(){
	return this->run_or_not;
}
