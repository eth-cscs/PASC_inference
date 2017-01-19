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

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_filename", boost::program_options::value< std::string >(), "name of input file with signal data (vector in PETSc format) [string]")
		("test_filename_out", boost::program_options::value< std::string >(), "name of output file with filtered signal data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of input file with original signal data without noise (vector in PETSc format) [string]")
		("test_filename_gamma0", boost::program_options::value< std::string >(), "name of input file with initial gamma approximation (vector in PETSc format) [string]")
		("test_save_all", boost::program_options::value<bool>(), "save results for all epssqr, not only for the best one [bool]")
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
		std::cout << "test_epssqr has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	int K, annealing; 
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not, save_all;

	std::string filename;
	std::string filename_out;
	std::string filename_solution;
	std::string filename_gamma0;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_filename", &filename, "data/samplesignal.bin");
	consoleArg.set_option_value("test_filename_out", &filename_out, "samplesignal");
	consoleArg.set_option_value("test_filename_solution", &filename_solution, "data/samplesignal_solution.bin");
	consoleArg.set_option_value("test_save_all", &save_all, false);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, false);
	consoleArg.set_option_value("test_printstats", &printstats, false);
	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, true);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, "shortinfo/samplesignal.txt");

	/* maybe gamma0 is given in console parameters */
	bool given_gamma0;
	if(consoleArg.set_option_value("test_filename_gamma0", &filename_gamma0)){
		given_gamma0 = true;
	} else {
		given_gamma0 = false;
	}

	/* maybe theta is given in console parameters */
	bool given_Theta;
	std::vector<double> Theta_list;
	double Theta_solution[K];
	if(consoleArg.set_option_value("test_Theta", &Theta_list)){
		given_Theta = true;
		
		/* control number of provided Theta */
		if(Theta_list.size() != K){
			coutMaster << "number of provided Theta solutions is different then number of clusters!\n";
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

	coutMaster << "- PROBLEM INFO ----------------------------\n";
	coutMaster << " DDT_size                    = " << std::setw(30) << DDT_size << " (decomposition in space)\n";
	coutMaster << " test_K                      = " << std::setw(30) << K << " (number of clusters)\n";
	if(given_Theta){
		coutMaster << " test_Theta                  = " << std::setw(30) << print_array(Theta_solution,K) << "\n";
	}

	coutMaster << " test_filename               = " << std::setw(30) << filename << " (name of input file with signal data)\n";
	coutMaster << " test_filename_out           = " << std::setw(30) << filename_out << " (name of output file with filtered signal data)\n";
	coutMaster << " test_filename_solution      = " << std::setw(30) << filename_solution << " (name of input file with original signal data without noise)\n";
	if(given_gamma0){
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << filename_gamma0 << " (name of input file with initial gamma approximation)\n";
	} else {
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << "NO" << " (name of input file with initial gamma approximation)\n";
	}
	coutMaster << " test_save_all               = " << std::setw(30) << save_all << " (save results for all epssqr, not only for the best one)\n";
	coutMaster << " test_epssqr                 = " << std::setw(30) << print_vector(epssqr_list) << " (penalty parameters)\n";
	coutMaster << " test_annealing              = " << std::setw(30) << annealing << " (number of annealing steps)\n";
	coutMaster << " test_cutgamma               = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})\n";
	coutMaster << " test_cutdata                = " << std::setw(30) << cutdata << " (cut data to {0,1})\n";
	coutMaster << " test_scaledata              = " << std::setw(30) << scaledata << " (scale data to {-1,1})\n";
	coutMaster << " test_printstats             = " << std::setw(30) << printstats << " (print basic statistics of data)\n";
	coutMaster << " test_shortinfo              = " << std::setw(30) << shortinfo_write_or_not << " (save shortinfo file after computation)\n";
	coutMaster << " test_shortinfo_header       = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)\n";
	coutMaster << " test_shortinfo_values       = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)\n";
	coutMaster << " test_shortinfo_filename     = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)\n";
	coutMaster << "-------------------------------------------\n" << "\n";


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
	coutMaster << "- start program\n";

/* 1.) prepare preliminary time-series data (to get the size of the problem T) */
	coutMaster << "--- PREPARING PRELIMINARY DATA ---\n";
	Signal1DData<PetscVector> mydata(filename);

/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---\n";

	/* prepare decomposition based on preloaded data */
	Decomposition decomposition(mydata.get_Tpreliminary(), 1, K, 1, DDT_size);

	/* print info about decomposition */
	decomposition.print(coutMaster);

/* 3.) prepare time-series data */
	coutMaster << "--- APPLY DECOMPOSITION TO DATA ---\n";
	mydata.set_decomposition(decomposition);

	/* print information about loaded data */
	mydata.print(coutMaster);

	/* print statistics */
	if(printstats) mydata.printstats(coutMaster);

/* 4.) prepare and load solution */
	Vec solution_Vec;
	TRYCXX( VecDuplicate(mydata.get_datavector()->get_vector(),&solution_Vec) );
	GeneralVector<PetscVector> solution(solution_Vec);
	solution.load_global(filename_solution);

/* 5.) prepare model */
	coutMaster << "--- PREPARING MODEL ---\n";

	/* prepare model on the top of given data */
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0]);

	/* print info about model */
	mymodel.print(coutMaster,coutAll);

/* 6.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---\n";

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* if gamma0 is provided, then load it */
	if(given_gamma0){
		coutMaster << " - loading and setting gamma0\n";
		mydata.load_gammavector(filename_gamma0);
	}

	/* print info about solver */
	mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
	/* write header of short output */
	if(shortinfo_write_or_not){
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;

		/* add info about the problem */
		oss_short_output_header << shortinfo_header << "K,epssqr,abserr,";

		/* append Theta solution */
		for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";

		/* print info from solver */
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		/* append end of line */
		oss_short_output_header << "\n";

		/* write to shortinfo file */
		shortinfo.write(oss_short_output_header.str());
			
		/* clear streams for next writing */
		oss_short_output_header.str("");
	}	
	
/* 6.) solve the problem with epssqrs and remember best solution */
	double epssqr, epssqr_best;
	double abserr; /* actual error */
	double abserr_best = std::numeric_limits<double>::max(); /* the error of best solution */

	Vec gammavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_gammavector()->get_vector(),&gammavector_best_Vec) );

	Vec thetavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_thetavector()->get_vector(),&thetavector_best_Vec) );
	
	/* go throught given list of epssqr */
	for(int depth = 0; depth < epssqr_list.size();depth++){
		epssqr = epssqr_list[depth];
		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr << " ---\n";

		/* set new epssqr */
		mymodel.set_epssqr(epssqr);

		/* cut data */
		if(cutdata) mydata.cutdata(0,1);

		/* scale data */
		if(scaledata){
			mydata.scaledata(-1,1,0,1);
		}
		
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before save */
		if(scaledata){
			mydata.scaledata(0,1,-1,1);
		}

		/* compute absolute error of computed solution */
		abserr = mydata.compute_abserr_reconstructed(solution);
		
		coutMaster << " - abserr = " << abserr << std::endl;
//		mysolver.printtimer(coutMaster);
//		mysolver.printstatus(coutMaster);	

		/* store obtained solution */
		if(save_all){
			coutMaster << "--- SAVING OUTPUT ---" << std::endl;
			oss << filename_out << "_epssqr" << epssqr;
			mydata.saveSignal1D(oss.str(),false);
			oss.str("");
		}
		
		/* store short info */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters and info about the problem */
			oss_short_output_values << shortinfo_values << K << "," << epssqr << "," << abserr << ",";

			/* append Theta solution */
			oss_short_output_values << mydata.print_thetavector(); 

			/* print info from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			oss_short_output_values << "\n";

			/* write to shortinfo file */
			shortinfo.write(oss_short_output_values.str());
			
			/* clear streams for next writing */
			oss_short_output_header.str("");
			oss_short_output_values.str("");
		}
	
		/* if this solution is better then previous, then store it */
		if(abserr < abserr_best){
			abserr_best = abserr;
			epssqr_best = epssqr;
			TRYCXX(VecCopy(mydata.get_gammavector()->get_vector(),gammavector_best_Vec));
			TRYCXX(VecCopy(mydata.get_thetavector()->get_vector(),thetavector_best_Vec));
		}
		
	}

	/* set best computed solution back to data */
	TRYCXX(VecCopy(gammavector_best_Vec,mydata.get_gammavector()->get_vector()));
	TRYCXX(VecCopy(thetavector_best_Vec, mydata.get_thetavector()->get_vector()));

/* 8.) store best solution */
	coutMaster << "--- SAVING OUTPUT ---" << std::endl;
	oss << filename_out;
	mydata.saveSignal1D(oss.str(),false);
	oss.str("");

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

