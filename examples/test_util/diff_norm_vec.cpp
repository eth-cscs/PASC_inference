/** @file diff_norm_vec.cpp
 *  @brief compute the difference of two PETSc vectors stored in given files
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

//typedef pascinference::algebra::PetscVector PetscVector;

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

	TRYCXX( VecCreate(PETSC_COMM_WORLD,&in1_Vec) );
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&in2_Vec) );

	GeneralVector<pascinference::algebra::PetscVector> in1(in1_Vec);
	GeneralVector<pascinference::algebra::PetscVector> in2(in2_Vec);

	in1.load_local(in1_filename);
	in2.load_local(in2_filename);

	/* print properties of vectors */
	coutMaster << std::setprecision(17);	
	coutMaster << std::endl;
	coutMaster << "in1:" << std::endl;
	coutMaster << " size      = " << std::setw(30) << in1.size() << std::endl;
	coutMaster << " norm      = " << std::setw(30) << norm(in1) << std::endl;
	coutMaster << " sum       = " << std::setw(30) << sum(in1) << std::endl;
	coutMaster << " max       = " << std::setw(30) << max(in1) << std::endl;
	coutMaster << " min       = " << std::setw(30) << min(in1) << std::endl;

	coutMaster << "in2:" << std::endl;
	coutMaster << " size      = " << std::setw(30) << in2.size() << std::endl;
	coutMaster << " norm      = " << std::setw(30) << norm(in2) << std::endl;
	coutMaster << " sum       = " << std::setw(30) << sum(in2) << std::endl;
	coutMaster << " max       = " << std::setw(30) << max(in2) << std::endl;
	coutMaster << " min       = " << std::setw(30) << min(in2) << std::endl;
	coutMaster << std::endl;


	double mynorm = norm(in1 - in2);
	coutMaster << " norm(in1 - in2) = " << mynorm << std::endl;

	Finalize();

	return 0;
}

