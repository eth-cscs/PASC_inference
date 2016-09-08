/** @file test_test_signal1D.cu
 *  @brief test the kmeans problem solver on simple 1D signal problem
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/signal1Ddata.h"
#include "model/graphh1fem.h"

#include <vector>

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
		("test_filename", boost::program_options::value< std::string >(), "name of input file with signal data (vector in PETSc format) [string]")
		("test_filename_out", boost::program_options::value< std::string >(), "name of output file with filtered signal data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of input file with original signal data without noise (vector in PETSc format) [string]")
		("test_epssqr", boost::program_options::value<std::vector<double> >()->multitoken(), "penalty parameters [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to set {0;1} [bool]")
		("test_scaledata", boost::program_options::value<bool>(), "scale to interval {-1,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to interval {0,1} [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_Theta", boost::program_options::value<std::vector<double> >()->multitoken(), "given solution Theta [K*int]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
		
	} else {
		std::cout << "test_epssqr has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	int K, annealing; 
	double graph_coeff; 
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not;

	std::string filename;
	std::string filename_out;
	std::string filename_solution;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_filename", &filename, "data/samplesignal.bin");
	consoleArg.set_option_value("test_filename_out", &filename_out, "samplesignal");
	consoleArg.set_option_value("test_filename_solution", &filename_solution, "");
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_printstats", &printstats, false);
	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, true);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, "shortinfo/samplesignal.txt");

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

	/* set decomposition in space */
	int DDT_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " DDT_size                    = " << std::setw(30) << DDT_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta                  = " << std::setw(30) << print_array(Theta_solution,K) << std::endl;
	}

	coutMaster << " test_filename               = " << std::setw(30) << filename << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_out           = " << std::setw(30) << filename_out << " (name of output file with filtered signal data)" << std::endl;
	coutMaster << " test_filename_solution      = " << std::setw(30) << filename_solution << " (name of input file with original signal data without noise)" << std::endl;
	coutMaster << " test_epssqr                 = " << std::setw(30) << print_vector(epssqr_list) << " (penalty parameters)" << std::endl;
	coutMaster << " test_annealing              = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma               = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata                = " << std::setw(30) << cutdata << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata              = " << std::setw(30) << scaledata << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_printstats             = " << std::setw(30) << printstats << " (print basic statistics of data)" << std::endl;
	coutMaster << " test_shortinfo              = " << std::setw(30) << shortinfo_write_or_not << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header       = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values       = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename     = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl << std::endl;


	/* start logging */
	std::ostringstream oss;
	oss << "log/" << filename_out << ".txt";
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

/* 1.) prepare preliminary time-series data (to get the size of the problem T) */
	coutMaster << "--- PREPARING PRELIMINARY DATA ---" << std::endl;
	Signal1DData<PetscVector> mydata(filename);

/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;

	/* prepare decomposition based on preloaded data */
	Decomposition decomposition(mydata.get_Tpreliminary(), 1, K, 1, DDT_size);

	/* print info about decomposition */
	decomposition.print(coutMaster);

/* 3.) prepare time-series data */
	coutMaster << "--- APPLY DECOMPOSITION TO DATA ---" << std::endl;

	mydata.set_decomposition(decomposition);

	/* print information about loaded data */
	mydata.print(coutMaster);

	/* print statistics */
	if(printstats) mydata.printstats(coutMaster);

	/* cut data */
	if(cutdata) mydata.cutdata(0,1);

	/* scale data */
	if(scaledata) mydata.scaledata(-1,1,0,1);

/* 4.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;

	/* prepare model on the top of given data */
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0]);

	/* print info about model */
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
	if(scaledata) mydata.scaledata(0,1,-1,1);

	coutMaster << "--- SAVING OUTPUT ---" << std::endl;
	oss << filename_out << "_depth0" << "_epssqr" << epssqr_list[0];
//	mydata.saveImage(oss.str(),true);
	oss.str("");

	/* write short output */
	if(shortinfo_write_or_not){

		/* add provided strings from console parameters and info about the problem */
		oss_short_output_header << shortinfo_header << "filename,K,depth,epssqr,";
		oss_short_output_values << shortinfo_values << filename << "," << K << ",0,0.0,"; 

		/* append Theta solution */
		for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";
		oss_short_output_values << mydata.print_thetavector(); 

		/* print info from solver */
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		/* append end of line */
		oss_short_output_header << std::endl;
		oss_short_output_values << std::endl;

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
		mysolver.set_annealing(1);

		/* scale data before computation */
		if(scaledata) mydata.scaledata(-1,1,0,1);

		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[depth] << " ---" << std::endl;
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before export */
		if(scaledata) mydata.scaledata(0,1,-1,1);

		coutMaster << "--- SAVING OUTPUT ---" << std::endl;
		oss << image_out << "_depth" << depth << "_epssqr" << epssqr_list[depth];
//		mydata.saveImage(oss.str(),false);
		oss.str("");
		
		/* write short output */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters */
			oss_short_output_values << shortinfo_values << filename << "," << K << "," << depth << "," << epssqr_list[depth] << ",";

			/* append Theta solution */
			oss_short_output_values << mydata.print_thetavector(); 

			/* append data from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			oss_short_output_values << std::endl;

			/* write data */
			shortinfo.write(oss_short_output_values.str());

			/* clear streams for next time */
			oss_short_output_values.str("");
		}
	}

	/* print solution */
	coutMaster << "--- THETA SOLUTION ---" << std::endl;
	mydata.print_thetavector(coutMaster);

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
	mysolver.printtimer(coutMaster);

	/* print short info */
	coutMaster << "--- FINAL SOLVER INFO ---" << std::endl;
	mysolver.printstatus(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

