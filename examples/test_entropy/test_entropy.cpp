/** @file test_entropy.cpp
 *  @brief test the new entropy model
 *
 *  @author Lukas Pospisil, Anna Marchenko
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#ifndef USE_DLIB
 #error 'This example is for DLIB'
#endif

#define DEFAULT_EPSSQR 1
#define DEFAULT_K 1
#define DEFAULT_KM 1
#define DEFAULT_FILENAME_IN "data/entropy_small_data.bin"
#define DEFAULT_FILENAME_OUT "entropy_small"
#define DEFAULT_SAVEALL false
#define DEFAULT_SAVERESULT true
#define DEFAULT_ANNEALING 1
#define DEFAULT_CUTGAMMA false
#define DEFAULT_SCALEDATA false
#define DEFAULT_CUTDATA false
#define DEFAULT_PRINTSTATS false
#define DEFAULT_PRINTINFO false
#define DEFAULT_SHORTINFO true
#define DEFAULT_SHORTINFO_FILENAME "shortinfo/entropy_small.txt"

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("ENTROPY EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_Km", boost::program_options::value<int>(), "number of moments [int]")
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with signal data (vector in PETSc format) [string]")
		("test_filename_out", boost::program_options::value< std::string >(), "name of output file with filtered signal data (vector in PETSc format) [string]")
		("test_filename_gamma0", boost::program_options::value< std::string >(), "name of input file with initial gamma approximation (vector in PETSc format) [string]")
		("test_saveall", boost::program_options::value<bool>(), "save results for all epssqr, not only for the best one [bool]")
		("test_saveresult", boost::program_options::value<bool>(), "save the solution [bool]")
		("test_epssqr", boost::program_options::value<std::vector<double> >()->multitoken(), "penalty parameters [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to set {0;1} [bool]")
		("test_scaledata", boost::program_options::value<bool>(), "scale to interval {-1,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to interval {0,1} [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]")
		("test_Theta", boost::program_options::value<std::vector<std::string> >()->multitoken(), "given solution Theta [K*int]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	/* start to measure computing time time */
    Timer timer_all;
    timer_all.start();

    /* load epssqr list */
	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
	} else {
		/* list is not given, add some value */
		epssqr_list.push_back(DEFAULT_EPSSQR);
	}

	int K, Km, annealing; 
	bool cutgamma, scaledata, cutdata, printstats, printinfo, shortinfo_write_or_not, saveall, saveresult;

	std::string filename_in;
	std::string filename_out;
	std::string filename_gamma0;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, DEFAULT_K);
	consoleArg.set_option_value("test_Km", &Km, DEFAULT_KM);
	consoleArg.set_option_value("test_filename_in", &filename_in, DEFAULT_FILENAME_IN);
	consoleArg.set_option_value("test_filename_out", &filename_out, DEFAULT_FILENAME_OUT);
	consoleArg.set_option_value("test_saveall", &saveall, DEFAULT_SAVEALL);
	consoleArg.set_option_value("test_saveresult", &saveresult, DEFAULT_SAVERESULT);
	consoleArg.set_option_value("test_annealing", &annealing, DEFAULT_ANNEALING);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, DEFAULT_CUTGAMMA);
	consoleArg.set_option_value("test_scaledata", &scaledata, DEFAULT_SCALEDATA);
	consoleArg.set_option_value("test_cutdata", &cutdata, DEFAULT_CUTDATA);
	consoleArg.set_option_value("test_printstats", &printstats, DEFAULT_PRINTSTATS);
	consoleArg.set_option_value("test_printinfo", &printinfo, DEFAULT_PRINTINFO);
	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, DEFAULT_SHORTINFO);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, DEFAULT_SHORTINFO_FILENAME);

	/* maybe gamma0 is given in console parameters */
	bool given_gamma0;
	if(consoleArg.set_option_value("test_filename_gamma0", &filename_gamma0)){
		given_gamma0 = true;
	} else {
		given_gamma0 = false;
	}

	/* maybe theta is given in console parameters */
	bool given_Theta;
	std::vector<std::string> Theta_list;
	double Theta_solution[K*Km];
	if(consoleArg.set_option_value("test_Theta", &Theta_list)){
		given_Theta = true;
		
		/* control number of provided Theta */
		if(Theta_list.size() != K){
			coutMaster << "number of provided Theta solutions is different then number of clusters! (you provided " << Theta_list.size() << " parameters)" << std::endl;
			return 0;
		}

		/* parse strings to doubles */
		if(!parse_strings_to_doubles(K,Km, Theta_list, Theta_solution) ){
			coutMaster << "unable to parse input Theta values!" << std::endl;
			return 0;
		}

	} else {
		given_Theta = false;
	}	

	/* set decomposition in space */
	int DDT_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
#ifdef USE_GPU
	coutMaster << " computing on GPU" << std::endl;
#else
	coutMaster << " computing on CPU" << std::endl;
#endif
	coutMaster << " DDT_size                    = " << std::setw(30) << DDT_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_Km                     = " << std::setw(30) << Km << " (number of moments)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta                  = " << std::setw(30) << print_array(Theta_solution,K,Km) << std::endl;
	}

	coutMaster << " test_filename_in            = " << std::setw(30) << filename_in << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_out           = " << std::setw(30) << filename_out << " (name of output file with filtered signal data)" << std::endl;
	if(given_gamma0){
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << filename_gamma0 << " (name of input file with initial gamma approximation)" << std::endl;
	} else {
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << "NO" << " (name of input file with initial gamma approximation)" << std::endl;
	}
	coutMaster << " test_saveall               = " << std::setw(30) << saveall << " (save results for all epssqr, not only for the best one)" << std::endl;
	coutMaster << " test_epssqr                 = " << std::setw(30) << print_vector(epssqr_list) << " (penalty parameters)" << std::endl;
	coutMaster << " test_annealing              = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma               = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata                = " << std::setw(30) << cutdata << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata              = " << std::setw(30) << scaledata << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_printstats             = " << std::setw(30) << printbool(printstats) << " (print basic statistics of data)" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(30) << printbool(printinfo) << " (print informations about created objects)" << std::endl;
	coutMaster << " test_shortinfo              = " << std::setw(30) << shortinfo_write_or_not << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header       = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values       = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename     = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;


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
	Signal1DData<PetscVector> mydata(filename_in);

/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;

	/* prepare decomposition based on preloaded data */
	Decomposition<PetscVector> decomposition(mydata.get_Tpreliminary(), 1, K, 1, DDT_size);

	/* print info about decomposition */
	if(printinfo) decomposition.print(coutMaster);

/* 3.) prepare time-series data */
	coutMaster << "--- APPLY DECOMPOSITION TO DATA ---" << std::endl;
	mydata.set_decomposition(decomposition);

	/* print information about loaded data */
	if(printinfo) mydata.print(coutMaster);

	/* print statistics */
	if(printstats) mydata.printstats(coutMaster);

/* 5.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;

	/* prepare model on the top of given data */
	EntropyH1FEMModel<PetscVector> mymodel(mydata, Km, epssqr_list[0]);

	/* print info about model */
	if(printinfo) mymodel.print(coutMaster,coutAll);

/* 6.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* if gamma0 is provided, then load it */
	if(given_gamma0){
		coutMaster << " - loading and setting gamma0" << std::endl;
		mydata.load_gammavector(filename_gamma0);
	}

	/* print info about solver */
	if(printinfo) mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
/* 6.) solve the problem with epssqrs and remember best solution */
	double epssqr;
	double epssqr_best = -1;
	double L; /* actual value of object function */
	double nbins;
	double L_best = std::numeric_limits<double>::max(); /* best solution */

	Vec gammavector_best_Vec; /* here we store solution with best L */
	TRYCXX( VecDuplicate(mydata.get_gammavector()->get_vector(),&gammavector_best_Vec) );

	Vec thetavector_best_Vec; /* here we store solution with best L */
	TRYCXX( VecDuplicate(mydata.get_thetavector()->get_vector(),&thetavector_best_Vec) );

	/* go throught given list of epssqr */
	for(int depth = 0; depth < epssqr_list.size();depth++){
		epssqr = epssqr_list[depth];
		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr << " ---" << std::endl;

		/* set new epssqr */
		mymodel.set_epssqr(epssqr);

		/* cut data */
		if(cutdata) mydata.cutdata(0,1);

		/* scale data */
		if(scaledata){
			mydata.scaledata(-1,1);
		}
		
		/* !!! solve the problem */
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before save */
		if(scaledata){
//			mydata.unscaledata(-1,1); //TODO: this is wrong
		}

		/* get function value */
		L = mysolver.get_L();
		nbins = mydata.compute_gammavector_nbins();

		coutMaster << " - L = " << L << ", nbins = " << nbins << std::endl;
//		mysolver.printtimer(coutMaster);
//		mysolver.printstatus(coutMaster);	

		/* store obtained solution */
		if(saveall){
			coutMaster << "--- SAVING OUTPUT ---" << std::endl;
			oss << filename_out << "_epssqr" << epssqr;
			mydata.saveSignal1D(oss.str(),false);
			oss.str("");
		}
		

		/* store short info */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters and info about the problem */
			if(depth==0) oss_short_output_header << shortinfo_header << "K,epssqr,L,";
			oss_short_output_values << shortinfo_values << K << "," << epssqr << "," << L << ",";
			
			/* append Theta solution */
			if(depth==0) for(int k=0; k < mymodel.get_thetavectorlength_local(); k++) oss_short_output_header << "Theta" << k << ",";
			oss_short_output_values << mydata.print_thetavector(); 

			/* append nbins */
			if(depth==0) oss_short_output_header << "Gamma nbins,";
			oss_short_output_values << nbins << ","; 

			/* print info from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			if(depth==0) oss_short_output_header << "\n";
			oss_short_output_values << "\n";

			/* write to shortinfo file */
			if(depth==0) shortinfo.write(oss_short_output_header.str());
			shortinfo.write(oss_short_output_values.str());
			
			/* clear streams for next writing */
			oss_short_output_header.str("");
			oss_short_output_values.str("");
		}
	
		/* if this solution is better then previous, then store it */
		if(L < L_best){
			L_best = L;
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
	coutMaster << " - with best epssqr = " << epssqr_best << std::endl;
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

	/* print info about power consumption */
	timer_all.stop();

	coutMaster << "--- FINAL INFO --------------------------------" << std::endl;
	coutMaster << "- total time: " << timer_all.get_value_last() << " s" << std::endl;
	coutMaster << "------------------------------------------------" << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

