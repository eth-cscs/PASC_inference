/** @file test_consoleoutput.cpp
 *  @brief test class and methods: ConsoleOutput 
 *
 *  Test the console output. Works without PETSC (one process) as well as with PETSC on several processes.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;

int main( int argc, char *argv[] ){

/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	if(!Initialize(argc, argv)){
		/* this happen for example when program is called with "--help" parameter */
		return 0;
	} 

/* --- TEST OUTPUT --- */
	/* say hello from master */
	coutMaster << "this is a test output from master (i.e. \"Hello World\" from master)" << std::endl;
	
	/* in this example, we use instance of GlobalManager, which could be used for manipulating with processes ids */
	coutMaster << "sample is running on " << GlobalManager.get_size() << " processes " << std::endl;

	/* say hello from this process */
	coutAll << "this a \"Hello World\" from pocess with id " << GlobalManager.get_rank() << std::endl;

	/* flush buffered output from all ranks to be sure that all of them finished cout */
	coutAll.synchronize();

	/* test output on more lines, the output has to be ended with std::endl */
	coutAll << "first line\nsecond line\nthird line" << std::endl;
	coutAll.synchronize();

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
