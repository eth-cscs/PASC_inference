/** @file util_print_vec.cpp
 *  @brief print the content of given PETSc vector
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_PRINT_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("in_filename", boost::program_options::value< std::string >(), "input vector [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string in_filename;

	if(!consoleArg.set_option_value("in_filename", &in_filename)){
		std::cout << "in_filename has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;;
	coutMaster << " in_filename            = " << std::setw(30) << in_filename << " (PETSc vector)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec in_Vec;

	TRYCXX( VecCreate(PETSC_COMM_WORLD,&in_Vec) );

	GeneralVector<PetscVector> in(in_Vec);

	in.load_local(in_filename);

	coutMaster << in << std::endl;

	Finalize<PetscVector>();

	return 0;
}

