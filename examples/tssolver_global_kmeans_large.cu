/** @file tssolver_global_kmeans_large.cu
 *  @brief test the varx global problem solver
 *
 *  Generate random Kmeans problem and solve it using VarX global solver.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver_global.h"
#include "data/tsdata_global.h"
#include "model/varxh1fem_global.h"

#include "kmeans3D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

extern int pascinference::DEBUG_MODE;

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

	/* read command line arguments */
	if(argc < 3){
		std::cout << "1. argument - T      - the dimension of the problem" << std::endl;
		std::cout << "2. argument - K      - number of clusters" << std::endl;
		std::cout << "3. argument - epssqr - penalty" << std::endl;
		std::cout << std::endl << argv[0] << " T K epssqr" << std::endl;
		return 1;
	}
	int k; /* iterators */
	int T = atoi(argv[1]); 
	int Kforall = atoi(argv[2]); 
	int epssqrforall = atof(argv[3]); 
	std::cout << "T       = " << T << " (length of time-series)" << std::endl;
	std::cout << "K       = " << Kforall << " (number of clusters)" << std::endl;
	std::cout << "epssqrt = " << epssqrforall << " (penalty)" << std::endl;
	
	Initialize(argc, argv); // TODO: load parameters from console input

	/* print info about environment */
	coutMaster << "---------------------" << std::endl;
	coutMaster << "nproc:   " << GlobalManager.get_size() << std::endl;
	coutAll <<    " my_rank: " << GlobalManager.get_rank() << std::endl;
	coutAll.synchronize();
	coutMaster << "---------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/kmeans_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int xdim = 3; /* data dimension */

	/* solution - for generating the problem */
	int solution_K = 3;
	int solution_xmem = 0;
	
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

	/* model parameters */
	int num = GlobalManager.get_size();
	int *K = new int(num);
	double *epssqr = new double(num);
	int *xmem = new int(num);

	int i;
	for(i=0; i<num; i++ ){
		K[i] = Kforall;
		epssqr[i] = epssqrforall;
		xmem[i] = 0;
	}
	
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel_Global mymodel(T, xdim, num, K, xmem, epssqr);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData_Global mydata(mymodel);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	mymodel.generate_data(solution_K, solution_xmem, solution_theta, solution_xstart, &solution_get_cluster_id, &mydata, false);
	mymodel.generate_data_add_noise(&mydata, noise_covariance, &solution_get_cluster_id);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver_Global mysolver(mydata);

	mysolver.setting.maxit = 1000;
	mysolver.setting.debug_mode = 2;
//	mysolver.print(coutMaster,coutAll);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans3D::saveVTK("results/kmeans",".vtk",T,num,K,mydata.get_datavector(),mydata.get_gammavector());
	coutAll.synchronize();
	
	/* print timers */
//	mysolver.printtimer(coutAll);
//	coutAll.synchronize();

	/* save results into CSV file */
	coutMaster << "--- SAVING CSV ---" << std::endl;
	mymodel.saveCSV("results/kmeans",&mydata);
	coutAll.synchronize();

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

