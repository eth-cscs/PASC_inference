/** @file edf.cu
 *  @brief test the varx global problem solver
 *
 *  Load EDF file and solve the problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/edfdata.h"
#include "model/graphh1fem.h"

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
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to {0,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to given interval [bool]")
		("test_cutdata_down", boost::program_options::value<double>(), "lower bound used for cutting data [double]")
		("test_cutdata_up", boost::program_options::value<double>(), "upper bound used for cutting data [double]")
		("test_max_record_nmb", boost::program_options::value<int>(), "maximum nuber of loaded records")
		("test_savevtk", boost::program_options::value<bool>(), "save results into vtk format [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K, max_record_nmb, annealing; 
	double epssqr; 
	bool cutgamma, savevtk, printstats, cutdata;
	double cutdata_up;
	double cutdata_down;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	consoleArg.set_option_value("test_max_record_nmb", &max_record_nmb, -1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_cutdata_down", &cutdata_down, -200);
	consoleArg.set_option_value("test_cutdata_up", &cutdata_up, 200);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_savevtk", &savevtk, true);
	consoleArg.set_option_value("test_printstats", &printstats, false);

	coutMaster << "----------------------- PROBLEM INFO --------------------------" << std::endl << std::endl;
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_epssqr             = " << std::setw(30) << epssqr << " (penalty)" << std::endl;
	coutMaster << " test_annealing          = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma           = " << std::setw(30) << cutgamma << " (cut gamma to {0,1})" << std::endl;
	coutMaster << " test_cutdata            = " << std::setw(30) << cutdata << " (cut data to given interval)" << std::endl;
	coutMaster << " test_cutdata_down       = " << std::setw(30) << cutdata_down << " (lower bound used for cutting data)" << std::endl;
	coutMaster << " test_cutdata_up         = " << std::setw(30) << cutdata_up << " (upper bound used for cutting data)" << std::endl;
	coutMaster << " test_max_record_nmb     = " << std::setw(30) << max_record_nmb << " (max number of loaded time-steps)" << std::endl;
	coutMaster << " test_savevtk            = " << std::setw(30) << savevtk << " (save results into vtk format)" << std::endl;
	coutMaster << " test_printstats         = " << std::setw(30) << printstats << " (print basic statistics of data)" << std::endl;
	coutMaster << "---------------------------------------------------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/edf_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	EdfData<PetscVector> mydata("data/S001R01.edf.txt", max_record_nmb);

	if(cutdata){
		mydata.cutdata(cutdata_down,cutdata_up);
	}

	if(printstats){
		mydata.printstats(coutMaster);
	}
	
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	BGMGraph mygraph("data/Koordinaten_EEG_P.bin");
	mygraph.process(2.5);
	GraphH1FEMModel<PetscVector> mymodel(mydata, mygraph, K, epssqr);

	mymodel.print(coutMaster,coutAll);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<PetscVector> mysolver(mydata, annealing);

	mysolver.print(coutMaster,coutAll);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma){
		mydata.cut_gamma();
	}

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

	/* print short info */
	coutMaster << "--- FINAL SOLVER INFO ---" << std::endl;
	mysolver.printstatus(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

