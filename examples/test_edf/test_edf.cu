/** @file edf.cu
 *  @brief test the varx global problem solver
 *
 *  Load EDF file and solve the problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/edfsolver.h"
#include "data/edfdata.h"
#include "model/edfh1fem.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

typedef petscvector::PetscVector PetscVector;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps")
		("test_max_record_nmb", boost::program_options::value<int>(), "maximum nuber of loaded records");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K, max_record_nmb, annealing; 
	double epssqr; 

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	consoleArg.set_option_value("test_max_record_nmb", &max_record_nmb, -1);
	consoleArg.set_option_value("test_annealing", &annealing, 1);

	coutMaster << "- PROBLEM INFO -----------------" << std::endl;
	coutMaster << " K              = " << K << " (number of clusters)" << std::endl;
	coutMaster << " epssqr         = " << epssqr << " (penalty)" << std::endl;
	coutMaster << " annealing      = " << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " max_record_nmb = " << max_record_nmb << " (max number of loaded time-steps)" << std::endl;
	coutMaster << "------------------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/edf_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	EdfData<PetscVector> mydata("data/S001R01.edf.txt", max_record_nmb);

//	mydata.print(coutMaster,coutAll);

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	BGM_Graph mygraph("data/Koordinaten_EEG_P.bin");
	mygraph.process(2.5);
	EdfH1FEMModel<PetscVector> mymodel(mydata, mygraph, K, epssqr);

	mymodel.print(coutMaster,coutAll);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	EdfSolver<PetscVector> mysolver(mydata, annealing);

	mysolver.print(coutMaster,coutAll);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* save results into CSV file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	mydata.saveVTK("edf_ascii");

	/* save results into CSV file */
//	coutMaster << "--- SAVING CSV ---" << std::endl;
//	mydata.saveCSV("results/edf.csv");

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
	mysolver.printtimer(coutAll);
	coutAll.synchronize();

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

