#include "pascinference.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#include "matrix/blockdiaglaplace_explicit.h"
#include "feasibleset/simplex.h"


using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
//typedef minlin::threx::HostVector<double> Global;

extern bool petscvector::PETSC_INITIALIZED;
extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
		
	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;
	
	/* say hello */	
	Message("- start program");

	/* dimension of the problem */
	int T = 5; /* length of time-series (size of the block) */
	int K = 2; /* number of clusters (block) */
	int N = K*T; /* global size */

/* ----------- SOLUTION IN PETSC -----------*/
	GeneralVector<Global> x(N); /* solution */

	GeneralVector<Global> x0(N); /* initial approximation */
	x0(gall) = 0.0;

	GeneralVector<Global> b(N); /* linear term */
	b(gall) = 1.0;
	
	BlockDiagLaplaceExplicitMatrix<Global> A(b,K); /* hessian matrix */

	SimplexFeasibleSet<Global> feasibleset(T,K); /* feasible set */

	/* add A,b,x0, to data */
	QPData<Global> data;
	data.set_A(&A);
	data.set_b(&b);
	data.set_x0(&x0);
	data.set_feasibleset(&feasibleset);
	data.set_x(&x);

	/* prepare solver */
	QPSolver<Global> myqp(data);

	/* solve the problem */
	pascinference::DEBUG_MODE = 10;
	myqp.solve(SOLVER_SPGQP);
	pascinference::DEBUG_MODE = 0;

	/* print solution */
	std::cout << *(data.get_x()) << std::endl;

	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

