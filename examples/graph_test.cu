/** @file graph_test.cu
 *  @brief test the graph matrix multiplication on several architectures
 *
 *  Generate n random vectors of length KN and compute matrix-multiplication.
 *  Measure the computing time
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "matrix/blockgraph.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
//		("test_filename", boost::program_options::value<char*>(), "name of file with coordinates")
		("test_T", boost::program_options::value<int>(), "dimension of the problem");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T;
//	char filename[50];
	
	consoleArg.set_option_value("test_T", &T, 10);
//	consoleArg.set_option_value("filename", &filename, "Koordinaten_EEG2.txt");

	coutMaster << "T        = " << std::setw(7) << T << " (length of time-series)" << std::endl;
//	coutMaster << "filename = " << std::setw(7) << filename << " (name of file with coordinates)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/projection_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;





	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

