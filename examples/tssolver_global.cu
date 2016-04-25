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
	int T = 11; /* length of time-series */

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
	int Knum = 2;
	int K[Knum] = {2,3}; /* number of clusters */
	Knum = GlobalManager.get_size(); // TODO: this is a hotfix for proc < Knum


	int xmem = 0; /* coeficient of Var model - memory of x */
	int umem = 0; /* coeficient of Var model - memory of u */
	
	double epssqr = 10; /* penalty */

	/* print info about environment */
	coutMaster << "---------------------" << std::endl;
	coutMaster << "nproc:   " << GlobalManager.get_size() << std::endl;
	coutAll <<    "my_rank: " << GlobalManager.get_rank() << std::endl;
	coutMaster << "---------------------" << std::endl << std::endl;
	
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel_Global mymodel(T, dimx, Knum, K, xmem, umem, epssqr);
	mymodel.print(coutMaster,coutAll);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData_Global mydata(mymodel);
	mydata.print(coutMaster);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::KMeans3D::generate(T, Solution_K, Solution_mu, Solution_diag_covariance, mydata.get_datavector());
	mydata.printcontent(coutMaster,coutAll);


	/* prepare time-series solver */
//	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
//	TSSolver<Global> mysolver(mydata);

	/* solve the problem */
//	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mysolver.setting.maxit = 10;
//	mysolver.setting.debug_mode = 10;
	
//	mysolver.print(coutMaster);

//	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans3D::saveVTK("output",".vtk",T,Knum,K,mydata.get_datavector(),mydata.get_gammavector());
	
//	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

