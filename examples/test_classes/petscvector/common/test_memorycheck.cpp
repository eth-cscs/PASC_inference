/** @file test_memorycheck.cpp
 *  @brief test class and methods: MemoryCheck
 *
 *  Test several utils for memory management. Could be used to get informations about the size of used memory.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;

/** the size of testing int array */
#define TEST_INTARRAY_SIZE 1000000

/* define testing function in separated namespace */
namespace test_namespace {
	/* print all known status of the memory */
	void print_memory_status(){
		coutMaster << "status (in bytes):" << std::endl;
		coutMaster << " Virtual memory all/used (percentage)  : ";
		coutMaster << MemoryCheck::get_virtual_all() << " / " << MemoryCheck::get_virtual_used() << " (" << MemoryCheck::get_virtual() << "%)" << std::endl;
		coutMaster << " Physical memory                       : ";
		coutMaster << MemoryCheck::get_physical() << std::endl;
	}

}

int main( int argc, char *argv[] ){

/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	if(!Initialize(argc, argv)){
		/* this happen for example when program is called with "--help" parameter */
		return 0;
	} 

/* --- TEST USED MEMORY --- */
	int old_value; /* here I store the physical memory */

	/* print status of the memory at the beginning of program */
	coutMaster << "Start program" << std::endl;
	coutMaster.push(); /* increase the offset of print */
	test_namespace::print_memory_status();
	coutMaster.pop(); /* decrease the offset of the print */
	coutMaster << std::endl;

	/* allocate some memory and see what happened */
	coutMaster << "Allocate array of " << TEST_INTARRAY_SIZE << " integers and fill it with sequence 0,1,2,..." << std::endl;
	old_value = MemoryCheck::get_virtual_used(); /* store old value of used physical memory */

	int *myarray = new int[TEST_INTARRAY_SIZE]; /* allocate array */
	for(int i=0;i<TEST_INTARRAY_SIZE;i++) myarray[i] = i; /* fill array */

	coutMaster.push();
	test_namespace::print_memory_status();
	coutMaster.pop();
	coutMaster << "It seems that the allocation took: " << MemoryCheck::get_virtual_used() - old_value << std::endl;
	coutMaster << std::endl;

	/* try to free the memory and see what changed in memory */
	coutMaster << "Free allocated array" << std::endl;
	free(myarray);
	coutMaster.push();
	test_namespace::print_memory_status();
	coutMaster.pop();
	coutMaster << "The difference before allocation and after free is: " << MemoryCheck::get_virtual_used() - old_value << std::endl;
	coutMaster << std::endl;

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
