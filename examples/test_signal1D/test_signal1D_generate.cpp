/** @file test_signal1D_generate.cpp
 *  @brief generate testing k-means sample
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

int myget_cluster_id(int t, int T){
	/* returns id of cluster based on t */
	double step = T/11.0;
	int cluster_id=0;
	if( ( (t >= 1*step) && (t < 3*step) ) || ( (t >= 4*step) && (t < 5*step) ) || ( (t >= 8*step) && (t < 9*step) ) ) {
		cluster_id = 1;
	}
	if( ( (t >= 3*step) && (t < 4*step) ) || ( (t >= 7*step) && (t < 8*step) ) || (t >= 9*step) ) {
		cluster_id = 2;
	}
	return cluster_id;
}


int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename", boost::program_options::value< std::string >(), "name of output file with signal data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of output file with original signal data without noise (vector in PETSc format) [string]")
		("test_filename_gamma0", boost::program_options::value< std::string >(), "name of output file with initial gamma approximation (vector in PETSc format) [string]")
		("test_T", boost::program_options::value< std::string >(), "length of 1D time-series [int]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	/* start to measure time */
    Timer timer_all;
    timer_all.start();

	int T;
	std::string filename;
	std::string filename_solution;
	std::string filename_gamma0;

	consoleArg.set_option_value("test_filename", &filename, "data/samplesignal.bin");
	consoleArg.set_option_value("test_filename_solution", &filename_solution, "data/samplesignal_solution.bin");
	consoleArg.set_option_value("test_filename_gamma0", &filename_solution, "data/samplesignal_gamma0.bin");
	consoleArg.set_option_value("test_T", &T, 10);

	/* set decomposition in space */
	int DDT_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " ranks_per_node              = " << std::setw(30) << ranks_per_node << " (number of MPI processes on one node)" << std::endl;
	coutMaster << " test_T                      = " << std::setw(30) << T << " (length of time-series)" << std::endl;
	coutMaster << " test_filename               = " << std::setw(30) << filename << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_solution      = " << std::setw(30) << filename_solution << " (name of input file with original signal data without noise)" << std::endl;
	coutMaster << " test_filename_gamma0        = " << std::setw(30) << filename_gamma0 << " (name of input file with initial gamma approximation)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;


	/* start logging */
	std::ostringstream oss;
	oss << "log/" << filename_out << ".txt";
	logging.begin(oss.str());
	oss.str("");

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* allocate vector of data */



	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

