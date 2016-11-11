/** @file stat_vec.cpp
 *  @brief print the properties of given PETSc vector
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

//typedef pascinference::algebra::PetscVector PetscVector;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_STAT_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("in_filename", boost::program_options::value< std::string >(), "input vector [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	std::string in_filename;

	if(!consoleArg.set_option_value("in_filename", &in1_filename)){
		std::cout << "in_filename has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	int K, annealing, width, height; 
	double graph_coeff; 
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not, graph_save;

	coutMaster << "- UTIL INFO ----------------------------\n";
	coutMaster << " in_filename            = " << std::setw(30) << in_filename << " (PETSc vector)\n";
	coutMaster << "-------------------------------------------\n" << "\n";

	/* prepare vector */
	Vec in_Vec;

	TRY( VecCreate(PETSC_COMM_WORLD,&in_Vec) );

	GeneralVector<pascinference::algebra::PetscVector> in(in_Vec);

	in.load_local(in_filename);

	/* print properties of vectors */
	coutMaster << std::setprecision(17);	
	coutMaster << std::endl;
	coutMaster << "properties of given vector:" << std::endl;
	coutMaster << " size      = " << std::setw(30) << in.size() << std::endl;
	coutMaster << " norm      = " << std::setw(30) << norm(in) << std::endl;
	coutMaster << " sum       = " << std::setw(30) << sum(in) << std::endl;
	coutMaster << " max       = " << std::setw(30) << max(in) << std::endl;
	coutMaster << " min       = " << std::setw(30) << min(in) << std::endl;

	Finalize();

	return 0;
}

