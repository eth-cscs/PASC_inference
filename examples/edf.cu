/** @file edf.cu
 *  @brief test the varx global problem solver
 *
 *  Load EDF file and solve the problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver_global.h"
#include "data/edfdata.h"
#include "model/varxh1fem_global.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K; 
	double epssqr; 

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);

	std::cout << "K       = " << K << " (number of clusters)" << std::endl;
	std::cout << "epssqrt = " << epssqr << " (penalty)" << std::endl;

	/* print info about environment */
	coutMaster << "---------------------" << std::endl;
	coutMaster << " nproc:   " << GlobalManager.get_size() << std::endl;
	coutAll <<    " my_rank: " << GlobalManager.get_rank() << std::endl;
	coutAll.synchronize();
	coutMaster << "---------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/edf_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int xdim = 3; /* data dimension */

	/* solution - for generating the problem */
	int solution_K = 3;
	int solution_xmem = 0;

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	EdfData mydata("data/S001R01.edf.txt");

	/* prepare model */
//	coutMaster << "--- PREPARING MODEL ---" << std::endl;
//	VarxH1FEMModel_Global mymodel(T, xdim, K, xmem, epssqr);

	/* prepare time-series solver */
//	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
//	TSSolver_Global mysolver(mydata);

//	mysolver.maxit = 1000;
//	mysolver.debug_mode = 2;
//	mysolver.print(coutMaster,coutAll);

	/* solve the problem */
//	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mysolver.solve();

	/* save results into CSV file */
//	coutMaster << "--- SAVING CSV ---" << std::endl;
//	mymodel.saveCSV("results/kmeans",&mydata);
//	coutAll.synchronize();

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
//	mysolver.printtimer(coutAll);
	coutAll.synchronize();

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

