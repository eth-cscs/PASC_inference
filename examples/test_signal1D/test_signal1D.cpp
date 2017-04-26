/** @file test_signal1D.cu
 *  @brief test the kmeans problem solver on simple 1D signal problem
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_fem_type", boost::program_options::value<int>(), "type of used FEM to reduce problem [0=FEM_SUM/1=FEM_HAT]")
		("test_fem_reduce", boost::program_options::value<double>(), "parameter of the reduction of FEM nodes [int,-1=false]")
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
		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]")
		("test_saveresult", boost::program_options::value<bool>(), "save the solution [bool]")
		("test_Theta", boost::program_options::value<std::vector<double> >()->multitoken(), "given solution Theta [K*int]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	/* start to measure power consumption and time */
    Timer timer_all;
    timer_all.start();
	const int ranks_per_node = PowerCheck::get_ranks_per_node();
    double node_energy    = PowerCheck::get_node_energy()/(double)ranks_per_node;
    double device_energy  = PowerCheck::get_device_energy()/(double)ranks_per_node;

    /* load epssqr list */
	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
		
	} else {
		std::cout << "test_epssqr has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	coutMaster << "kunda?" << std::endl;


	int K, annealing, fem_type; 
	bool cutgamma, scaledata, cutdata, printstats, printinfo, shortinfo_write_or_not, save_all, saveresult;
	double fem_reduce;

	std::string filename;
	std::string filename_out;
	std::string filename_solution;
	std::string filename_gamma0;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_fem_type", &fem_type, 1);
	consoleArg.set_option_value("test_fem_reduce", &fem_reduce, 1.0);
	consoleArg.set_option_value("test_filename", &filename, "data/samplesignal.bin");
	consoleArg.set_option_value("test_filename_out", &filename_out, "samplesignal");
	consoleArg.set_option_value("test_filename_solution", &filename_solution, "data/samplesignal_solution.bin");
	consoleArg.set_option_value("test_save_all", &save_all, false);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, false);
	consoleArg.set_option_value("test_printstats", &printstats, false);
	consoleArg.set_option_value("test_printinfo", &printinfo, false);
	consoleArg.set_option_value("test_saveresult", &saveresult, true);
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
#ifdef USE_CUDA
	coutMaster << " computing on GPU" << std::endl;
#else
	coutMaster << " computing on CPU" << std::endl;
#endif
	coutMaster << " ranks_per_node              = " << std::setw(30) << ranks_per_node << " (number of MPI processes on one node)" << std::endl;
	coutMaster << " DDT_size                    = " << std::setw(30) << DDT_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta                  = " << std::setw(30) << print_array(Theta_solution,K) << std::endl;
	}

	coutMaster << " test_fem_type               = " << std::setw(30) << fem_type << " (type of used FEM to reduce problem [0=FEM_SUM/1=FEM_HAT])" << std::endl;
	coutMaster << " test_fem_reduce             = " << std::setw(30) << fem_reduce << " (parameter of the reduction of FEM node)" << std::endl;
	coutMaster << " test_filename               = " << std::setw(30) << filename << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_out           = " << std::setw(30) << filename_out << " (name of output file with filtered signal data)" << std::endl;
	coutMaster << " test_filename_solution      = " << std::setw(30) << filename_solution << " (name of input file with original signal data without noise)" << std::endl;
	if(given_gamma0){
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << filename_gamma0 << " (name of input file with initial gamma approximation)" << std::endl;
	} else {
		coutMaster << " test_filename_gamma0        = " << std::setw(30) << "NO" << " (name of input file with initial gamma approximation)" << std::endl;
	}
	coutMaster << " test_save_all               = " << std::setw(30) << save_all << " (save results for all epssqr, not only for the best one)" << std::endl;
	coutMaster << " test_epssqr                 = " << std::setw(30) << print_vector(epssqr_list) << " (penalty parameters)" << std::endl;
	coutMaster << " test_annealing              = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma               = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata                = " << std::setw(30) << cutdata << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata              = " << std::setw(30) << scaledata << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_saveresult             = " << std::setw(30) << saveresult << " (save reconstructed signal)" << std::endl;
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
	Signal1DData<PetscVector> mydata(filename);

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

/* 4.) prepare and load solution */
	Vec solution_Vec;
	TRYCXX( VecDuplicate(mydata.get_datavector()->get_vector(),&solution_Vec) );
	GeneralVector<PetscVector> solution(solution_Vec);
	solution.load_global(filename_solution);

/* 5.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;

	/* prepare FEM reduction */
	Fem<PetscVector> *fem;
	if(fem_type == 0){
		fem = new Fem<PetscVector>(fem_reduce);
	}
	if(fem_type == 1){
		fem = new FemHat<PetscVector>(fem_reduce);
	}

	/* prepare model on the top of given data */
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0], fem);

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
	double abserr; /* actual error */
	double abserr_best = std::numeric_limits<double>::max(); /* the error of best solution */

	Vec gammavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_gammavector()->get_vector(),&gammavector_best_Vec) );

	Vec thetavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_thetavector()->get_vector(),&thetavector_best_Vec) );

	/* energy for one iteration */
	double node_energy_it;
    	double node_energy_it_sum;
	
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
			mydata.scaledata(-1,1,0,1);
		}
		
		/* measure energy at begin */
		MPI_Barrier(MPI_COMM_WORLD);
		node_energy_it    = PowerCheck::get_node_energy()/(double)ranks_per_node;

		/* !!! solve the problem */
		mysolver.solve();

		/* measure energy in the end */
		MPI_Barrier(MPI_COMM_WORLD);
		node_energy_it     = PowerCheck::get_node_energy()/(double)ranks_per_node - node_energy_it;
		node_energy_it_sum = PowerCheck::mpi_sum_reduce(node_energy_it);
		
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
		if(save_all && saveresult){
			coutMaster << "--- SAVING OUTPUT ---" << std::endl;
			oss << filename_out << "_epssqr" << epssqr;
			mydata.saveSignal1D(oss.str(),false);
			oss.str("");
		}
		

		/* store short info */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters and info about the problem */
			if(depth==0) oss_short_output_header << shortinfo_header << "K,epssqr,abserr,energy,";
			oss_short_output_values << shortinfo_values << K << "," << epssqr << "," << abserr << "," << node_energy_it_sum << ",";
			
			/* append Theta solution */
			if(depth==0) for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";
			oss_short_output_values << mydata.print_thetavector(); 

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
	if(saveresult){
		coutMaster << "--- SAVING OUTPUT ---" << std::endl;
		coutMaster << " - with best epssqr = " << epssqr_best << std::endl;
		oss << filename_out;
		mydata.saveSignal1D(oss.str(),false);
		oss.str("");
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

	/* print info about power consumption */
	timer_all.stop();
	MPI_Barrier(MPI_COMM_WORLD);
	node_energy    = PowerCheck::get_node_energy()/(double)ranks_per_node - node_energy;
	device_energy  = PowerCheck::get_device_energy()/(double)ranks_per_node - device_energy;

	double total_node_energy = PowerCheck::mpi_sum_reduce(node_energy);
	double total_device_energy = PowerCheck::mpi_sum_reduce(device_energy);
	double time_all = timer_all.get_value_sum();

	coutMaster << "--- ENERGY INFO --------------------------------" << std::endl;
	coutMaster << "- total time: " << time_all << " s" << std::endl;
	#ifdef USE_CRAYPOWER
		coutMaster << "- node" << std::endl;
		coutMaster.push();
		coutMaster << "- " << total_node_energy << " Joules, ";
		coutMaster << total_node_energy/time_all << " Watts";
		coutMaster << std::endl;
		coutMaster.pop();
		coutMaster << "- device" << std::endl;
		coutMaster.push();
		coutMaster << "- " << total_device_energy << " Joules, ";
		coutMaster << total_device_energy/time_all << " Watts";
		coutMaster << std::endl;
		coutMaster.pop();
		coutMaster << "- node + device" << std::endl;
		coutMaster.push();
		coutMaster << "- " << (total_node_energy + total_device_energy) << " Joules, ";
		coutMaster << (total_node_energy + total_device_energy)/time_all << " Watts";
		coutMaster << std::endl;
		coutMaster.pop();
	#endif

	coutMaster << "------------------------------------------------" << std::endl;



	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

