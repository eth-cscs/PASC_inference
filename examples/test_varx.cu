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

#include "varx2D.h"

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

	/* dimension of the problem */
	int dimx = 2; /* data dimension */
	int T = 1000; /* length of time-series */
	int K = 2; /* number of clusters */
	int xmem = 3; /* coeficient of Var model */
	int umem = 1; /* coeficient of Var model */

	double epssqr = 10; /* penalty */
	
/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Global> mymodel(T, dimx, K, xmem, umem, epssqr);
//	mymodel.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);
//	mydata.print(coutMaster);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::Varx2D<Global>::generate(T,xmem,mydata.get_datavector(),mydata.get_u());
//	mydata.printcontent(coutMaster);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

	mymodel.get_gammadata()->print(coutMaster);
	mymodel.get_thetadata()->printcontent(coutMaster);


//	mysolver.setting.maxit = 50;
//	mysolver.setting.debug_mode = 0;

	/* solve the problem */
	/* gamma_solver = SOLVER_SPGQP, theta_solver = SOLVER_CG */
//	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mysolver.solve(SOLVER_SPGQP, SOLVER_CG);

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::Varx2D<Global>::saveVTK("output.vtk",T+xmem,K,mydata.get_datavector(),mydata.get_gammavector());
	
//	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

