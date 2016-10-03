/** @file test_offset.cpp
 *  @brief test class and methods: Offset
 *
 *  Test manipulation with offset in coutMaster and coutAll.
 *  Please note that coutMaster and coutAll uses global instance of Offset called offset.
 *  Therefore they use shared offset and push()/pop() is applied to both of them simultaneously.
 *  (i.e. it doesn't matter if you call coutMaster.pop() or coutAll.pop(), in both cases the global offset in increased )
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;


int main( int argc, char *argv[] ){

/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	if(!Initialize(argc, argv)){
		return 0;
	} 

/* --- TEST OFFSET --- */
	/* print something pretty */
	coutMaster << "I have to go into shop and buy:" << std::endl;
	coutMaster.push();
	coutMaster << "- something to drink" << std::endl;
	coutMaster.push();
	coutMaster << "- coke, sprite, soda" << std::endl;
	coutMaster << "- something alcoholic" << std::endl;
	coutMaster.push();
	coutMaster << "- rum" << std::endl << "- vodka" << std::endl << "- beer" << std::endl;
	coutMaster.pop();
	coutMaster << "- non-alcoholic wine" << std::endl;
	coutMaster.pop();
	coutMaster << "- something to eat" << std::endl;
	coutMaster.push();
	coutMaster << "- grilled lama" << std::endl << "- baked crocodile" << std::endl;
	coutMaster.pop();
	coutMaster.pop();
	coutMaster << std::endl;

	coutMaster << "So, which process will go to shop?" << std::endl;
	coutMaster.push();
	/* go through processes and ask them. If I am the process, then I will answer */
	for(int i=0; i < GlobalManager.get_size();i++){ 
		coutMaster << "What about you, process number " << i << "?" << std::endl;
	
		if(i == GlobalManager.get_rank()){
			/* I should say something */
			coutAll.push();
			coutAll << "me not, because:" << std::endl;
			coutAll.push();
			coutAll << "I am lazy" << std::endl << "I don't have a time" << std::endl << "it is too late" << std::endl;
			coutAll.pop();
			coutAll.pop();
		}
	}
	coutAll.synchronize(); /* has to be called to flush print from all processes */
	coutMaster.pop();

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
