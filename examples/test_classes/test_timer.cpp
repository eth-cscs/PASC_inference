/** @file test_timer.cpp
 *  @brief test class and methods: Timer and StackTimer
 *
 *  Test time measurement.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

/* the size of problem */
#define TEST_UNIT_SIZE 1e4

using namespace pascinference;

/* define testing functions in separated namespace */
namespace test_namespace {

/* the function which takes some time */
int test_dosomething(int number=1000) {
	LOG_FUNC_STATIC_BEGIN

	/* exponential complexity */
	long int sum;
	for(int i=0;i<number;i++){
		sum = 0;
		for(int j=0;j<i;j++){
			sum += j;
		}
	}
	
	LOG_FUNC_STATIC_END
	
	return sum;
}

} /* end of namespace */

int main( int argc, char *argv[] ){

/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	if(!Initialize(argc, argv)){
		/* this happen for example when program is called with "--help" parameter */
		return 0;
	} 

/* --- TEST TIMERS --- */
	coutMaster << "- testing timers" << std::endl;
	
	/* define variables */
	Timer mytimer, mytimer2;
	StackTimer mystacktimer;

	/* initialize Timer */
	mytimer.restart();
	mytimer2.restart();

	coutMaster << "        size         |      Timer (sub)     |        Timer         |    StackTimer (sub)  |       StackTimer      " << std::endl;
	double mystacktimer_value, mystacktimer_value2; /* to stop StackTimer in "exactly" same time as Timer */
	for(int i=TEST_UNIT_SIZE;i<=10*TEST_UNIT_SIZE;i+=TEST_UNIT_SIZE){
		/* start to measure time */
		mytimer.start();
		mystacktimer.start();

		/* perform some action to measure time */
		test_namespace::test_dosomething(i);

			/* while timers are running, perform the subaction and measure also this subaction */
			mytimer2.start(); /* start the timer of subaction */
			mystacktimer.start(); /* the advantage of StackTimer - I don't need a new timer, I just increase the depth of stack */
			test_namespace::test_dosomething(TEST_UNIT_SIZE);
			mytimer2.stop(); /* stop the timer of subaction */
			mystacktimer_value2 = mystacktimer.stop(); /* it removes the timer from stack */

		/* stop to measure time of outer iteration */
		mytimer.stop();
		mystacktimer_value = mystacktimer.stop(); /* it removes the timer from stack */

		/* std::setw(int) is used to set fixed length of following variable output */
		coutMaster << std::setw(20) << i << " | ";


		/* get the value of time between Timer.start() and Timer.stop() */
		coutMaster << std::setw(20) << mytimer2.get_value_last() << " | ";
		coutMaster << std::setw(20) << mytimer.get_value_last() << " | ";
		
		/* stacktimer works differently, when I call StackTimer.stop(), it returns the value and remove time from stack */
		coutMaster << std::setw(20) << mystacktimer_value2 << " | "; 
		coutMaster << std::setw(20) << mystacktimer_value << std::endl; 
	}

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
