/** @file tssolver_global.cu
 *  @brief test the varx global problem solver
 *
 *  Generate random kmeans 3D problem and solve it using VarX global solver.
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

int main( int argc, char *argv[] )
{
	
	Initialize(argc, argv); // TODO: load parameters from console input
	
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int dimx = 3; /* data dimension */
	int T = 10000; /* length of time-series */

	/* solution - for generating the problem */
	int Solution_K = 4;
	double Solution_mu[12] = {0.0, 0.0, 0.0,
					1.0, 0.0, 0.0,
					0.0, 1.0, 0.0,
					0.5, 0.5, 0.2};

	double Solution_diag_covariance[12] = {	0.05, 0.05, 0.01,  
									0.05, 0.05, 0.01,  
									0.05, 0.05, 0.01,  
									0.05,  0.05, 0.05};


	/* model parameters */
	int num = 5;
//	int K[Knum] = {4,3,1,2,5}; /* number of clusters */
	int K[num] = {1,2,3,4,5}; 
	double epssqr[num] = {20,20,20,20,20}; /* penalty */

	num = GlobalManager.get_size(); // TODO: this is a hotfix for proc < num


	int xmem = 0; /* coeficient of Var model - memory of x */
	int umem = 0; /* coeficient of Var model - memory of u */
	

	/* print info about environment */
	coutMaster << "---------------------" << std::endl;
	coutMaster << "nproc:   " << GlobalManager.get_size() << std::endl;
	coutAll <<    " my_rank: " << GlobalManager.get_rank() << std::endl;
	coutAll.synchronize();
	coutMaster << "---------------------" << std::endl << std::endl;
	
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel_Global mymodel(T, dimx, num, K, xmem, umem, epssqr);
	mymodel.print(coutMaster,coutAll);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData_Global mydata(mymodel);
	mydata.print(coutMaster);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::KMeans3D::generate(T, Solution_K, Solution_mu, Solution_diag_covariance, mydata.get_datavector());

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver_Global mysolver(mydata);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.setting.maxit = 20;
	mysolver.setting.debug_mode = 5;
	mysolver.print(coutMaster,coutAll);

	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans3D::saveVTK("output",".vtk",T,num,K,mydata.get_datavector(),mydata.get_gammavector());
	
	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

