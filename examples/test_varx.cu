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

	/* test - for solution of gamma problem */
//	GeneralVector<Global> gamma_solution = *(mydata.get_gammavector());
//	GeneralVector<Global> theta_solution = *(mydata.get_thetavector());

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::Varx2D<Global>::generate(T,xmem,0.00005, mydata.get_datavector(),mydata.get_u());
//	example::Varx2D<Global>::generate(T,xmem,0.00005, mydata.get_datavector(),mydata.get_u(),&theta_solution,&gamma_solution);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

//	mymodel.get_gammadata()->print(coutMaster);
//	mymodel.get_thetadata()->print(coutMaster);
//	mydata.printcontent(coutMaster);
	
	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.setting.maxit = 10;
	mysolver.setting.debug_mode = 10;
	
	/* test: give solution of gamma and see what will happen */
//	*(mydata.get_gammavector()) = gamma_solution;
//	*(mydata.get_thetavector()) = theta_solution;

	mysolver.solve();

	mydata.printcontent(coutMaster);


	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::Varx2D<Global>::saveVTK("output.vtk",T,xmem,K,mydata.get_datavector(),mydata.get_gammavector());
	
//	mysolver.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

