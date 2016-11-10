/** @file diff_norm_vec.cpp
 *  @brief compute the difference of two PETSc vectors stored in given files
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

typedef petscvector::PetscVector PetscVector;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_DIFF_NORM_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("in1_filename", boost::program_options::value< std::string >(), "first input vector [string]")
		("in2_filename", boost::program_options::value< std::string >(), "second input vector [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	std::string in1_filename;
	std::string in2_filename;

	if(!consoleArg.set_option_value("in1_filename", &in1_filename) || !consoleArg.set_option_value("in2_filename", &in2_filename)){
		std::cout << "in1_filename and in2_filename have to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	int K, annealing, width, height; 
	double graph_coeff; 
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not, graph_save;

	coutMaster << "- UTIL INFO ----------------------------\n";
	coutMaster << " in1_filename            = " << std::setw(30) << in1_filename << " (first vector)\n";
	coutMaster << " in2_filename            = " << std::setw(30) << in2_filename << " (second vector)\n";
	coutMaster << "-------------------------------------------\n" << "\n";

	/* prepare vectors */
	Vec in1_Vec;
	Vec in2_Vec;

	GeneralVector<PetscVector> in1(in1_Vec);
	GeneralVector<PetscVector> in2(in2_Vec);

	in1.load_global(in1_filename);
	in2.load_global(in2_filename);

	double mynorm = norm(in1 - in2);
	coutMaster << std::setprecision(17);	
	coutMaster << " norm(in1 - in2) = " << mynorm << std::endl;

	Finalize();

	return 0;
}

