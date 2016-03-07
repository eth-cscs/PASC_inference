#include "pascinference.h"
#include "matrix/laplace_explicit_regular.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"
#include "result/qpresult.h"

using namespace pascinference;

/* set what is what ( which type of vector to use where) */
#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector Global;
	extern bool petscvector::PETSC_INITIALIZED;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> Host;
	typedef minlin::threx::DeviceVector<double> Device;
#endif

extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{

	
	/* say hello */	
	Message("- start program");

	/* dimension of the problem */
	int N = 10;

/* ----------- SOLUTION IN PETSC -----------*/
#ifdef USE_PETSCVECTOR
	Message("------- Solution in Petsc -----------");

	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;

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

	/* print some funny info */
	std::cout << data << std::endl;
	std::cout << result << std::endl;
	std::cout << myqp << std::endl;
	std::cout << x << std::endl;

	petscvector::PETSC_INITIALIZED = false;
	Finalize();
#endif

/* ----------- SOLUTION IN MINLIN -----------*/
#ifdef USE_MINLIN
	Message("------- Solution in Minlin -----------");


	GeneralVector<Host> xh(N); /* solution */

	GeneralVector<Host> x0h(N); /* initial approximation */
	xh(gall) = 0.0;

	GeneralVector<Host> bh(N); /* linear term */
	bh(gall) = 1.0;
	
	LaplaceExplicitRegularMatrix<Host> Ah(bh); /* hessian matrix */

	/* add A,b,x0 to data */
	QPData<Host> datah;
	datah.A = &Ah;
	datah.b = &bh;
	datah.x0 = &x0h;

	/* add x to results */
	QPResult<Host> resulth;
	resulth.x = &xh;

	QPSolver<Host> myqph(datah,resulth);

	myqph.solve(SOLVER_CG);

	/* print some funny info */
	std::cout << datah << std::endl;
	std::cout << resulth << std::endl;
	std::cout << myqph << std::endl;
	std::cout << xh << std::endl;

#endif

	/* say bye */	
	Message("- end program");

	return 0;
}

