/** @file benchmark2.cu
 *  @brief test the varx problem solver
 *
 *  Generate random kmeans 3D problem and solve it using VarX solver.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "model/varxh1fem.h"

#include "kmeans3D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
//typedef minlin::threx::HostVector<double> Host;

extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
	
	Initialize(argc, argv); // TODO: load parameters from console input
	
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int dimx = 3; /* data dimension */
	int T = 1000; /* length of time-series */
	int K = 4; /* number of clusters */
	int xmem = 0; /* coeficient of Var model - memory of x */
	int umem = 0; /* coeficient of Var model - memory of u */

	double mu[12] = {0.0, 0.0, 0.0,
					1.0, 0.0, 0.0,
					0.0, 1.0, 0.0,
					0.5, 0.5, 0.2};

	double diag_covariance[12] = {	0.05, 0.05, 0.01,  
									0.05, 0.05, 0.01,  
									0.05, 0.05, 0.01,  
									0.05,  0.05, 0.05};
	
	double epssqr = 10; /* penalty */
	
/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Global> mymodel(T, dimx, K, xmem, umem, epssqr);
	mymodel.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::KMeans3D<Global>::generate(T, 4, mu, diag_covariance, mydata.get_datavector());

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.setting.maxit = 10;
	mysolver.setting.debug_mode = 10;
	
	mysolver.print(coutMaster);

	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans3D<Global>::saveVTK("output.vtk",T,K,mydata.get_datavector(),mydata.get_gammavector());
	
	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

