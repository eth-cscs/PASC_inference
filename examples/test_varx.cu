/** @file test_varx.cu
 *  @brief test the varx problem solver
 *
 *  Generate random problem and solve it.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "model/varxh1fem.h"

#include "varx1D.h"

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
	int dim = 1; /* data dimension */
	int T = 30; /* length of time-series */
	int K = 3; /* number of clusters */

	/* parameters of the model */
	double mu[3] = {2.0, 1.0, 1.5};

	double dev[3] = {0.1,0.02,0.15};
	
	int Q = 1;
	double penalty = 10;
	
/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Global> mymodel(T, dim, K, Q, penalty);
	coutMaster << "Model: " << mymodel << std::endl;

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);
	coutMaster << "Data:  " << mydata << std::endl;

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::Varx1D<Global>::generate(T,K,mu,dev,mydata.get_datavector());
	mydata.printcontent(coutMaster);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

//	mysolver.setting.maxit = 50;
//	mysolver.setting.debug_mode = 0;

	/* solve the problem */
	/* gamma_solver = SOLVER_SPGQP, theta_solver = SOLVER_CG */
//	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mysolver.solve(SOLVER_SPGQP, SOLVER_CG);

	/* save results into VTK file */
//	coutMaster << "--- SAVING VTK ---" << std::endl;
//	example::KMeans2D<Global>::saveVTK("output.vtk",T,K,mydata.get_datavector(),mydata.get_gammavector());

//	mysolver.printtimer(coutMaster);

	/* say bye */	
//	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

