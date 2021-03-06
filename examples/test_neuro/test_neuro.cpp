/** @file test_neuro.cpp
 *  @brief test the neuro filtering possibilities
 *
 *  Load EDF file and solve the filtering problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

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
		("test_graph_save", boost::program_options::value<bool>(), "save VTK with graph or not [bool]")
		("test_data_out", boost::program_options::value< std::string >(), "part of output filename [string]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_Theta", boost::program_options::value<std::vector<double> >()->multitoken(), "given solution Theta [K*int]")
		("test_epssqr", boost::program_options::value<std::vector<double> >()->multitoken(), "penalty parameters [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_savevtk", boost::program_options::value<bool>(), "save results into vtk format [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to {0,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to given interval [bool]")
		("test_cutdata_down", boost::program_options::value<double>(), "lower bound used for cutting data [double]")
		("test_cutdata_up", boost::program_options::value<double>(), "upper bound used for cutting data [double]")
		("test_scaledata", boost::program_options::value<bool>(), "scale data to interval 0,1 [bool]")
		("test_shiftdata", boost::program_options::value<bool>(), "shift data with -0.5 [bool]")
		("test_shiftdata_coeff", boost::program_options::value<double>(), "coeficient of data shift [double]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
		
	} else {
		std::cout << "test_epssqr has to be set! Call application with parameter --h to see all parameters" << std::endl;
		return 0;
	}

	int K, max_record_nmb, annealing, DDT_size, DDR_size; 
	bool cutgamma, savevtk, printstats, cutdata, scaledata, shiftdata, shortinfo_write_or_not, graph_save;
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
	consoleArg.set_option_value("test_graph_save", &graph_save, false);
	
	consoleArg.set_option_value("test_data_out", &data_out, "test_edf");

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_savevtk", &savevtk, false);
	consoleArg.set_option_value("test_printstats", &printstats, false);

	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_cutdata_down", &cutdata_down, -200);
	consoleArg.set_option_value("test_cutdata_up", &cutdata_up, 200);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_shiftdata", &shiftdata, false);
	consoleArg.set_option_value("test_shiftdata_coeff", &shiftdata_coeff, 200);

	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, true);
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


	coutMaster << "----------------------------------- PROBLEM INFO --------------------------------------" << std::endl << "" << std::endl;
	coutMaster << " nmb of proc             = " << std::setw(30) << GlobalManager.get_size() << " (number of MPI processes)" << std::endl;
	coutMaster << " test_DDT                = " << std::setw(30) << DDT_size << " (decomposition in time)" << std::endl;
	coutMaster << " test_DDR                = " << std::setw(30) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << "" << std::endl;
	coutMaster << " test_data_filename      = " << std::setw(30) << data_filename << " (name of input file)" << std::endl;
	coutMaster << " test_max_record_nmb     = " << std::setw(30) << max_record_nmb << " (max number of loaded time-steps)" << std::endl;
	coutMaster << " test_graph_coordinates  = " << std::setw(30) << graph_coordinates << " (name of input file with coordinates)" << std::endl;
	coutMaster << " test_graph_coeff        = " << std::setw(30) << graph_coeff << " (threshold coefficient of graph)" << std::endl;
	coutMaster << " test_graph_save         = " << std::setw(30) << graph_save << " (save VTK with graph or not)" << std::endl;
	coutMaster << " test_data_out           = " << std::setw(30) << data_out << " (part of output filename)" << std::endl;
	coutMaster << "" << std::endl;
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta              = " << std::setw(30) << print_array(Theta_solution,K) << "" << std::endl;
	}
	coutMaster << " test_epssqr             = " << std::setw(30) << print_vector(epssqr_list) << " (penalty)" << std::endl;
	coutMaster << " test_annealing          = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_savevtk            = " << std::setw(30) << savevtk << " (save results into vtk format)" << std::endl;
	coutMaster << " test_printstats         = " << std::setw(30) << printstats << " (print basic statistics of data)" << std::endl;
	coutMaster << "" << std::endl;
	coutMaster << " test_cutgamma           = " << std::setw(30) << cutgamma << " (cut gamma to {0,1})" << std::endl;
	coutMaster << " test_cutdata            = " << std::setw(30) << cutdata << " (cut data to given interval)" << std::endl;
	coutMaster << " test_cutdata_down       = " << std::setw(30) << cutdata_down << " (lower bound used for cutting data)" << std::endl;
	coutMaster << " test_cutdata_up         = " << std::setw(30) << cutdata_up << " (upper bound used for cutting data)" << std::endl;
	coutMaster << " test_scaledata          = " << std::setw(30) << scaledata << " (scale data to interval [0,1])" << std::endl;
	coutMaster << " test_shiftdata          = " << std::setw(30) << shiftdata << " (shift data by coeficient)" << std::endl;
	coutMaster << " test_shiftdata_coeff    = " << std::setw(30) << shiftdata_coeff << " (shifting coeficient)" << std::endl;
	coutMaster << "" << std::endl;
	coutMaster << " test_shortinfo          = " << std::setw(30) << shortinfo_write_or_not << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header   = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values   = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "---------------------------------------------------------------------------------------" << std::endl << "" << std::endl;

	/* control the decomposition */
	if(DDT_size*DDR_size != GlobalManager.get_size()){
		coutMaster << "Sorry, DDT*DDR != nproc" << std::endl;
		return 0;
	}

	/* start logging */
	std::ostringstream oss;
	oss << "log/" << data_out << "_r" << max_record_nmb << "_K" << K << "_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss.str());
	oss.str("");
		
	/* start shortinfo output */
	if(shortinfo_write_or_not){
		shortinfo.begin(shortinfo_filename);
	}
	std::ostringstream oss_short_output_values;
	std::ostringstream oss_short_output_header;
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

/* 1.) prepare time-series data */
	coutMaster << "--- PREPARING PRELIMINARY DATA ---" << std::endl;
	EdfData<PetscVector> mydata(data_filename, max_record_nmb);

/* 2a.) prepare graph */
	coutMaster << "--- PREPARING GRAPH ---" << std::endl;
	BGMGraph<PetscVector> graph(graph_coordinates);
	graph.process(graph_coeff);
	graph.print(coutMaster);
	if(graph_save){
		/* save decoposed graph to see if space (graph) decomposition is working */
		oss << "results/" << data_out << "_DDR" << DDR_size << "_graph.vtk";
		graph.saveVTK(oss.str());
		oss.str("");
	}


/* 2b.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;
	Decomposition<PetscVector> mydecomposition(mydata.get_Tpreliminary(), graph, K, 1, DDT_size, DDR_size);
	mydecomposition.print(coutMaster);

/* 2c.) set new decomposition to data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	mydata.set_decomposition(mydecomposition);
	mydata.print(coutMaster);

	/* cut data to given interval [cutdata_down,cutdata_up] */
	if(cutdata)	mydata.cutdata(cutdata_down,cutdata_up);

	/* print basic stats about loaded data */
	if(printstats) mydata.printstats(coutMaster);
	
	/* scale data to interval [0,1] */
	if(scaledata) mydata.scaledata(0,1,cutdata_down,cutdata_up);

/* 3.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0]);
	mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* print info about solver */
	mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
/* 6.) solve the problem with initial epssqr */
	coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[0] << " ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma) mydata.cutgamma();

	/* unscale data before save */
	if(scaledata) mydata.unscaledata(0,1);


	coutMaster << "--- SAVING OUTPUT ---" << std::endl;
	oss << data_out << "_epssqr" << epssqr_list[0];
	mydata.saveVector(oss.str(),true);
	oss.str("");


	/* write short output */
	if(shortinfo_write_or_not){
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;
		oss_short_output_values << shortinfo_values;

		/* add info about the problem */
		oss_short_output_header << "max_record_nmb,K,depth,epssqr,";
		oss_short_output_values << max_record_nmb << "," << K << ",0," << epssqr_list[0] << ","; 

		/* append Theta solution */
		for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";
		oss_short_output_values << mydata.print_thetavector(); 

		/* print info from solver */
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		/* append end of line */
		oss_short_output_header << "" << std::endl;
		oss_short_output_values << "" << std::endl;

		/* write to shortinfo file */
		shortinfo.write(oss_short_output_header.str());
		shortinfo.write(oss_short_output_values.str());
			
		/* clear streams for next writing */
		oss_short_output_header.str("");
		oss_short_output_values.str("");
	
	}


/* 7.) solve the problems with other epssqr */
	for(int depth = 1; depth < epssqr_list.size();depth++){
		/* set new epssqr */
		mymodel.set_epssqr(epssqr_list[depth]);

		/* decrease the number of annealing steps in TSSolver to 1 */
//		mysolver.set_annealing(1);

		/* scale data before computation */
		if(scaledata) mydata.scaledata(0,1,cutdata_down,cutdata_up);

		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[depth] << " ---" << std::endl;
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before export */
		if(scaledata) mydata.unscaledata(0,1);

		coutMaster << "--- SAVING OUTPUT ---" << std::endl;
		oss << data_out << "_epssqr" << epssqr_list[depth];
		mydata.saveVector(oss.str(),false);
		oss.str("");

		/* write short output */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters */
			oss_short_output_values << shortinfo_values << max_record_nmb << "," << K << "," << depth << "," << epssqr_list[depth] << ",";

			/* append Theta solution */
			oss_short_output_values << mydata.print_thetavector(); 

			/* append data from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			oss_short_output_values << "" << std::endl;

			/* write data */
			shortinfo.write(oss_short_output_values.str());

			/* clear streams for next time */
			oss_short_output_values.str("");
		}
	}

/* 6.) save VTK - there is possibility to save only one problem into VTK */ 
	//TODO: make it in different way
	if(savevtk) {
		coutMaster << "--- SAVING VTK ---" << std::endl;
		mydata.saveVTK(data_out);
	}

/* 7.) print results */

	/* print solution */
	coutMaster << "--- THETA SOLUTION ---" << std::endl;
	mydata.print_thetavector(coutMaster);

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
	Finalize<PetscVector>();

	return 0;
}

