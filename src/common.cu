#include "common.h"

/* instance of global timer */
Timer timer;

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


}

/*!
 * final call of the application
 */ 
void Finalize(){

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





/* ------------ TIMER ------------ */

double Timer::getUnixTime(void){
	struct timespec tv;
	if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
	return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}

void Timer::start(){
	double t_start = getUnixTime();
	this->time_stack.push(t_start);
}
	
double Timer::stop(){
	double t_end = getUnixTime();
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
