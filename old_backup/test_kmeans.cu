/** @file test_kmeans.cu
 *  @brief test the kmeans problem solver
 *
 *  Generate random problem and solve it.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "model/kmeansh1fem.h"

#include "kmeans2D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
//typedef minlin::threx::HostVector<double> Global;

extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
	
	Initialize(argc, argv); // TODO: load parameters from console input
	
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* dimension of the problem */
	int dim = 2; /* data dimension */
	int T = 200; /* length of time-series */
	int K = 3; /* number of clusters */

	/* parameters of the model */
	double muK1[2] = {0.25, 0};
	double muK2[2] = {0.0, -0.5};
	double muK3[2] = {0.0, 0.5};
	double *mu[3] = {muK1,muK2,muK3};

	double covarianceK1[4] = {0.001, 0.0, 0.0, 0.1};
	double covarianceK2[4] = {0.005, 0.0, 0.0, 0.05};
	double covarianceK3[4] = {0.005, 0.0, 0.0, 0.05};
	double *covariance[3] = {covarianceK1,covarianceK2,covarianceK3};
	
	double penalty = 10;
	
/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	KmeansH1FEMModel<Global> mymodel(T, dim, K, penalty);
	coutMaster << "Model: " << mymodel << std::endl;

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);
	coutMaster << "Data:  " << mydata << std::endl;

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::KMeans2D<Global>::generate(T,K,mu,covariance,mydata.get_datavector());
//	mydata.printcontent(coutMaster);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

	mysolver.setting.maxit = 50;
	mysolver.setting.debug_mode = 0;

	/* solve the problem */
	/* gamma_solver = SOLVER_SPGQP, theta_solver = SOLVER_CG */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve(SOLVER_SPGQP, SOLVER_CG);

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans2D<Global>::saveVTK("output.vtk",T,K,mydata.get_datavector(),mydata.get_gammavector());

	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

