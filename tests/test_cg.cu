#include "pascinference.h"
#include "matrix/laplace_explicit_regular.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"
#include "result/qpresult.h"

using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
typedef minlin::threx::HostVector<double> Host;
typedef minlin::threx::DeviceVector<double> Device;

extern bool petscvector::PETSC_INITIALIZED;
extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
		
	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;
	
	/* say hello */	
	Message("- start program");

	/* dimension of the problem */
	int N = 10;

/* ----------- SOLUTION IN PETSC -----------*/
	GeneralVector<Global> x(N); /* solution */

	GeneralVector<Global> x0(N); /* initial approximation */
	x0(gall) = 0.0;

	GeneralVector<Global> b(N); /* linear term */
	b(gall) = 1.0;
	
	LaplaceExplicitRegularMatrix<Global> A(b); /* hessian matrix */

	/* add A,b,x0 to data */
	QPData<Global> data;
	data.A = &A;
	data.b = &b;
	data.x0 = &x0;

	/* add x to results */
	QPResult<Global> result;
	result.x = &x;

	QPSolver<Global> myqp(data,result);

	myqp.solve(SOLVER_CG);

	std::cout << *(result.x) << std::endl;

	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

