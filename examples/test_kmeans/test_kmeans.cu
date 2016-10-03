/** @file test_kmeans_large.cu
 *  @brief test the varx global problem solver with Kmeans problem
 *
 *  Generate random Kmeans problem and solve it using VarX global solver.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/kmeansdata.h"
#include "model/kmeansh1fem.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

int solution_get_cluster_id(int t, int T){
    int id_cluster;
    if((t >= ceil(0*T) && t < ceil(0.1*T)) || (t >= ceil(0.35*T) && t < ceil(0.4*T)) || (t >= ceil(0.5*T) && t < ceil(0.52*T)) || (t >= ceil(0.8*T) && t < ceil(0.9*T)) || (t >= ceil(0.2*T) && t < ceil(0.25*T))){
        id_cluster = 0;
    }

    if((t >= ceil(0.1*T) && t < ceil(0.2*T)) || (t >= ceil(0.4*T) && t < ceil(0.5*T)) || (t >= ceil(0.8*T) && t < ceil(0.8*T)) || (t >= ceil(0.95*T) && t <= ceil(1.0*T)) || (t >= ceil(0.6*T) && t < ceil(0.8*T))){
	id_cluster = 1;
    }
    
    if((t >= ceil(0.25*T) && t < ceil(0.35*T)) || (t >= ceil(0.52*T) && t < ceil(0.6*T)) || (t >= ceil(0.9*T) && t < ceil(0.95*T))){
	id_cluster = 2;
    }
    return id_cluster;
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int T, K, annealing;
	double epssqr; 

	consoleArg.set_option_value("test_T", &T, 1000);
	consoleArg.set_option_value("test_K", &K, 3);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	consoleArg.set_option_value("test_annealing", &annealing, 1);

	coutMaster << "- PROBLEM INFO -----------------\n";
	coutMaster << " T         = " << T << " (length of time-series)\n";
	coutMaster << " K         = " << K << " (number of clusters)\n";
	coutMaster << " epssqr    = " << epssqr << " (penalty)\n";
	coutMaster << " annealing = " << annealing << " (number of annealing steps)\n";
	coutMaster << "------------------------------\n" << "\n";

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/kmeans_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program\n";

	/* parameters of the model */
	int xdim = 3; /* data dimension */

	/* solution - for generating the problem */
	int solution_K = 3;
	
	double solution_theta[9] = {
		 0.0, 0.0, 0.0,		/* K=1,n=1,2,3:  mu */
		 1.0, 0.0, 0.0,		/* K=2,n=1,2,3 */
		 0.0, 1.0, 0.0		/* K=3,n=1,2,3 */
	};

	double solution_xstart[3] = {
		0.0, 0.0, 0.0	/* n=1,2,3 */
	};

	double noise_covariance[9] = {
		0.05, 0.05, 0.01,  
		0.05, 0.05, 0.01,  
		0.05, 0.05, 0.01  
	};

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---\n";
	KmeansData<PetscVector> mydata(T,xdim);

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---\n";
	KmeansH1FEMModel<PetscVector> mymodel(mydata, xdim, K, epssqr);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---\n";
	mydata.generate(solution_K, solution_theta, &solution_get_cluster_id, false);
	mydata.add_noise(noise_covariance);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---\n";
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---\n";
	mysolver.solve();

	/* print solution */
	coutMaster << "--- THETA SOLUTION ---\n";
	mydata.print_thetavector(coutMaster);

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---\n";
	mydata.saveVTK("results/test_kmeans.vtk");

	/* save results into CSV file */
	coutMaster << "--- SAVING CSV ---\n";
	mydata.saveCSV("results/test_kmeans.csv");

	/* print timers */
	coutMaster << "--- TIMERS INFO ---\n";
	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program\n";

	logging.end();
	Finalize();

	return 0;
}

