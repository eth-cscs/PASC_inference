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
		("filename_in", boost::program_options::value< std::string >(), "input vector [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string filename_in;

	if(!consoleArg.set_option_value("filename_in", &filename_in)){
		std::cout << "filename_in has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;;
	coutMaster << " filename_in            = " << std::setw(30) << filename_in << " (PETSc vector)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec in_Vec;

	TRYCXX( VecCreate(PETSC_COMM_WORLD,&in_Vec) );

	GeneralVector<PetscVector> in(in_Vec);

	in.load_local(filename_in);

	coutMaster << in << std::endl;

	Finalize<PetscVector>();

	return 0;
}

