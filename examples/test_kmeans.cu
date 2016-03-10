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
#include "model/kmeansh1.h"

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
	Message("- start program");

	/* dimension of the problem */
	int dim = 2; /* data dimension */
	int T = 5; /* length of time-series (size of the block) */
	int K = 3; /* number of clusters (block) */

/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	KmeansH1Model<Global> mymodel(T, dim, K);

	/* prepare time-series data */
	TSData<Global> mydata(mymodel);

	/* prepare time-series solver */
	TSSolver<Global> mysolver(mydata);
	std::cout << mysolver << std::endl;

	/* solve the problem */
	/* gamma_solver = SOLVER_SPGQP, theta_solver = SOLVER_CG */
	mysolver.solve(SOLVER_SPGQP, SOLVER_CG);

	/* say bye */	
	Message("- end program");
	
	Finalize();

	return 0;
}

