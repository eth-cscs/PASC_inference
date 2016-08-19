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
		("test_DDT", boost::program_options::value<int>(), "decomposition in time [int]")
		("test_DDR", boost::program_options::value<int>(), "decomposition in space [int]")
		("test_data_filename", boost::program_options::value< std::string >(), "name of input file [string]")
		("test_max_record_nmb", boost::program_options::value<int>(), "maximum nuber of loaded records")
		("test_graph_coordinates", boost::program_options::value< std::string >(), "name of input file with coordinates [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold coefficient of graph [double]")
		("test_data_out", boost::program_options::value< std::string >(), "part of output filename [string]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_Theta", boost::program_options::value<std::vector<double> >()->multitoken(), "given solution Theta [K*int]")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_savevtk", boost::program_options::value<bool>(), "save results into vtk format [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to {0,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to given interval [bool]")
		("test_cutdata_down", boost::program_options::value<double>(), "lower bound used for cutting data [double]")
		("test_cutdata_up", boost::program_options::value<double>(), "upper bound used for cutting data [double]")
		("test_scaledata", boost::program_options::value<bool>(), "scale data to interval -1,1 [bool]")
		("test_shiftdata", boost::program_options::value<bool>(), "shift data with -0.5 [bool]")
		("test_shiftdata_coeff", boost::program_options::value<double>(), "coeficient of data shift [double]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K, max_record_nmb, annealing, DDT_size, DDR_size; 
	double epssqr; 
	bool cutgamma, savevtk, printstats, cutdata, scaledata, shiftdata, shortinfo;
	double cutdata_up, cutdata_down, shiftdata_coeff, graph_coeff;

	std::string data_filename;
	std::string graph_coordinates;
	std::string data_out;

	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_DDT", &DDT_size, GlobalManager.get_size());
	consoleArg.set_option_value("test_DDR", &DDR_size, 1);

	consoleArg.set_option_value("test_data_filename", &data_filename, "data/S001R01.edf");
	consoleArg.set_option_value("test_max_record_nmb", &max_record_nmb, -1);
	consoleArg.set_option_value("test_graph_coordinates", &graph_coordinates, "data/Koordinaten_EEG_P.bin");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 2.5);
	consoleArg.set_option_value("test_data_out", &data_out, "test_edf");

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_epssqr", &epssqr, 0.000001);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_savevtk", &savevtk, true);
	consoleArg.set_option_value("test_printstats", &printstats, false);

	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_cutdata_down", &cutdata_down, -200);
	consoleArg.set_option_value("test_cutdata_up", &cutdata_up, 200);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_shiftdata", &shiftdata, false);
	consoleArg.set_option_value("test_shiftdata_coeff", &shiftdata_coeff, 200);

	consoleArg.set_option_value("test_shortinfo", &shortinfo, true);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, "shortinfo/myshortinfo.txt");

	/* maybe theta is given in console parameters */
	bool given_Theta;
	std::vector<double> Theta_list;
	double Theta_solution[K];
	if(consoleArg.set_option_value("test_Theta", &Theta_list)){
		given_Theta = true;
		
		/* control number of provided Theta */
		if(Theta_list.size() != K){
			coutMaster << "number of provided Theta solutions is different then number of clusters!" << std::endl;
			return 0;
		}

		/* store solution in array */
		for(int k=0;k < K;k++){
			Theta_solution[k] = Theta_list[k];
		}
	} else {
		given_Theta = false;
	}	


	coutMaster << "----------------------------------- PROBLEM INFO --------------------------------------" << std::endl << std::endl;
	coutMaster << " nmb of proc             = " << std::setw(30) << GlobalManager.get_size() << " (number of MPI processes)" << std::endl;
	coutMaster << " test_DDT                = " << std::setw(30) << DDT_size << " (decomposition in time)" << std::endl;
	coutMaster << " test_DDR                = " << std::setw(30) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << std::endl;
	coutMaster << " test_data_filename      = " << std::setw(30) << data_filename << " (name of input file)" << std::endl;
	coutMaster << " test_max_record_nmb     = " << std::setw(30) << max_record_nmb << " (max number of loaded time-steps)" << std::endl;
	coutMaster << " test_graph_coordinates  = " << std::setw(30) << graph_coordinates << " (name of input file with coordinates)" << std::endl;
	coutMaster << " test_graph_coeff        = " << std::setw(30) << graph_coeff << " (threshold coefficient of graph)" << std::endl;
	coutMaster << " test_data_out           = " << std::setw(30) << data_out << " (part of output filename)" << std::endl;
	coutMaster << std::endl;
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta              = " << std::setw(30) << "[";
		for(int k=0;k<K;k++){
			coutMaster << Theta_solution[k];
			if(k<K-1){
				coutMaster << ",";
			}
		}
		coutMaster << "]" << std::endl;
	}
	coutMaster << " test_epssqr             = " << std::setw(30) << epssqr << " (penalty)" << std::endl;
	coutMaster << " test_annealing          = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_savevtk            = " << std::setw(30) << savevtk << " (save results into vtk format)" << std::endl;
	coutMaster << " test_printstats         = " << std::setw(30) << printstats << " (print basic statistics of data)" << std::endl;
	coutMaster << std::endl;
	coutMaster << " test_cutgamma           = " << std::setw(30) << cutgamma << " (cut gamma to {0,1})" << std::endl;
	coutMaster << " test_cutdata            = " << std::setw(30) << cutdata << " (cut data to given interval)" << std::endl;
	coutMaster << " test_cutdata_down       = " << std::setw(30) << cutdata_down << " (lower bound used for cutting data)" << std::endl;
	coutMaster << " test_cutdata_up         = " << std::setw(30) << cutdata_up << " (upper bound used for cutting data)" << std::endl;
	coutMaster << " test_scaledata          = " << std::setw(30) << scaledata << " (scale data to interval [-1,1])" << std::endl;
	coutMaster << " test_shiftdata          = " << std::setw(30) << shiftdata << " (shift data by coeficient)" << std::endl;
	coutMaster << " test_shiftdata_coeff    = " << std::setw(30) << shiftdata_coeff << " (shifting coeficient)" << std::endl;
	coutMaster << std::endl;
	coutMaster << " test_shortinfo          = " << std::setw(30) << shortinfo << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header   = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values   = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "---------------------------------------------------------------------------------------" << std::endl << std::endl;

	/* control the decomposition */
	if(DDT_size*DDR_size != GlobalManager.get_size()){
		coutMaster << "Sorry, DDT*DDR != nproc" << std::endl;
		return 0;
	}

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "log/" << data_out << "_r" << max_record_nmb << "_K" << K << "_epssqr" << epssqr << "_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

/* 1.) prepare time-series data */
	coutMaster << "--- PREPARING PRELIMINARY DATA ---" << std::endl;
	EdfData<PetscVector> mydata(data_filename, max_record_nmb);

/* 2a.) prepare graph */
	coutMaster << "--- PREPARING GRAPH ---" << std::endl;
	BGMGraph mygraph(graph_coordinates);
	mygraph.process(graph_coeff);
	mygraph.print(coutMaster);

/* 2b.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;
	Decomposition mydecomposition(mydata.get_Tpreliminary(), mygraph, K, 1, DDT_size, DDR_size);
	mydecomposition.print(coutMaster);

/* 2c.) set new decomposition to data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	mydata.set_decomposition(mydecomposition);
	mydata.print(coutMaster);

	/* cut data to given interval [cutdata_down,cutdata_up] */
	if(cutdata){
		mydata.cutdata(cutdata_down,cutdata_up);
	}

	/* print basic stats about loaded data */
	if(printstats){
		mydata.printstats(coutMaster);
	}
	
	/* scale data to interval [-1,1] */
	if(scaledata){
		mydata.scaledata(-1,1);
	}

/* 3.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr);
	mymodel.print(coutMaster,coutAll);

/* 4.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<PetscVector> mysolver(mydata, annealing);
	mysolver.print(coutMaster,coutAll);

/* 5.) solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma){
		mydata.cutgamma();
	}

	/* unscale data to original interval */
	if(scaledata){
		mydata.unscaledata(-1,1);
	}

/* 6.) save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	mydata.saveVTK(data_out);

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

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

