/** @file test_globalmanager.cpp
 *  @brief test class and methods: GlobalManagerClass and GlobalManager
 *
 *  Test the GlobalManagerClass and its global instance GlobalManager. 
 *  This stuff can be used for manipulation with MPI_Comm_rank and MPI_Comm_size.
 *  Works also without PETSC - in that case size=1 and rank=0.
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
	/* get the size of COMM (i.e. number of processes) on Master and afterwards on each process */
	coutMaster << "master   GlobalManager.get_size() = " << GlobalManager.get_size() << ", GlobalManager.get_rank() = " << GlobalManager.get_rank() << std::endl;
	coutAll    << "all(me)  GlobalManager.get_size() = " << GlobalManager.get_size() << ", GlobalManager.get_rank() = " << GlobalManager.get_rank() << std::endl;
	coutAll.synchronize();

	/* say hello from all processes 0,1,... */
	coutMaster << std::endl;
	coutMaster << "testing FOR cycle with rank comparision" << std::endl;
	for(int i=0; i < GlobalManager.get_size();i++){
		/* if i is equal to my rank, than I will print */
		if(i==GlobalManager.get_rank()){
			coutAll << "i=" << i << ": I am printing because for me: GlobalManager.get_rank()=" << GlobalManager.get_rank() << std::endl;
		}
	}
	coutAll.synchronize();

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
