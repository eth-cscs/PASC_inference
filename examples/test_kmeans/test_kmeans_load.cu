/** @file test_kmeans_load.cu
 *  @brief 
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include "solver/kmeanssolver.h"
#include "data/kmeansdata.h"
#include "model/kmeansh1fem.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_data_filename", boost::program_options::value< std::string >(), "name of file with data (vector in PETSc format) [string]")
		("test_gamma0_filename", boost::program_options::value< std::string >(), "name of file with gamma0 (vector in PETSc format) [string]")
		("test_xdim", boost::program_options::value<int>(), "data dimension [int]")
		("test_K", boost::program_options::value<int>(), "number of clusters to test [int]")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter [double]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* get and print console parameters */
	std::string data_filename;
	std::string gamma0_filename;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;
	bool shortinfo;
	int K, xdim; 
	double epssqr;

	consoleArg.set_option_value("test_data_filename", &data_filename, "data/data_kmeans_T100.bin");
	consoleArg.set_option_value("test_gamma0_filename", &gamma0_filename, "data/gamma0_kmeans_T100K3.bin");
	consoleArg.set_option_value("test_xdim", &xdim, 3);
	consoleArg.set_option_value("test_K", &K, 3);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	consoleArg.set_option_value("test_shortinfo", &shortinfo, true);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, "results/kmeans_shortinfo.txt");


	coutMaster << "----------------------- PROBLEM INFO --------------------------" << std::endl << std::endl;
	coutMaster << " test_data_filename      = " << std::setw(30) << data_filename << " (name of file with data)" << std::endl;
	coutMaster << " test_gamma0_filename    = " << std::setw(30) << gamma0_filename << " (name of file with gamma0)" << std::endl;
	coutMaster << " test_xdim               = " << std::setw(30) << xdim << " (data dimension)" << std::endl;
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters to test)" << std::endl;
	coutMaster << " test_epssqp             = " << std::setw(30) << epssqr << " (penalty parameter)" << std::endl;
	coutMaster << " test_shortinfo          = " << std::setw(30) << shortinfo << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header   = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values   = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	
	coutMaster << "---------------------------------------------------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_filename_log;
	oss_filename_log << "results/test_kmeans_load_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_filename_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* prepare time-series data */
	coutMaster << "--- LOADING DATA ---" << std::endl;
	KmeansData<PetscVector> mydata(data_filename, xdim);

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	KmeansH1FEMModel<PetscVector> mymodel(mydata, xdim, K, epssqr);

	mydata.print(coutMaster,coutAll);
	
	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	KmeansSolver<PetscVector> mysolver(mydata);

	/* load gammavector from file */
	mydata.load_gammavector(gamma0_filename);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	mydata.saveVTK("results/test_kmeans.vtk");

	/* save results into CSV file */
	coutMaster << "--- SAVING CSV ---" << std::endl;
	mydata.saveCSV("results/test_kmeans.csv");

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* write short output */
	if(shortinfo){
		std::ostringstream oss_short_output_values;
		std::ostringstream oss_short_output_header;
		
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;
		oss_short_output_values << shortinfo_values;
		
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		std::ofstream myfile;
		myfile.open(shortinfo_filename.c_str());
		
		/* master writes the file with short info (used in batch script for quick computation) */
		if( GlobalManager.get_rank() == 0){
			myfile << oss_short_output_header.str() << std::endl;
			myfile << oss_short_output_values.str() << std::endl;
		}
		TRY( PetscBarrier(NULL) );
		
		myfile.close();
	}


	logging.end();
	Finalize();

	return 0;
}

