#include "pascinference.h"

#include "matrix/laplace_explicit.h"
#include "problem/qpproblem.h"


using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
typedef minlin::threx::HostVector<double> Host;
typedef minlin::threx::DeviceVector<double> Device;

extern bool petscvector::PETSC_INITIALIZED;
extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
	int N = 5; /* dimension of the problem */	
		
	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;
	
	/* say hello */	
	Message("- start program");

	/* prepare solution vector */
	GeneralVector<Global> x_global(N);
	x_global(gall) = 0.0;

	/* prepare initial approximation */
	GeneralVector<Global> x0_global(x_global);
	x0_global(gall) = 0.0;

	/* prepare right hand-side vector */
	GeneralVector<Global> b_global(x_global);
	b_global(gall) = 1.0/N;

	/* prepare matrix */	
	LaplaceExplicitMatrix<Global> A_global(x_global);

	/* prepare problem and fill it */
	QPProblem<Global> myqp;
	
	
//pascinference::DEBUG_MODE = 100;	

	std::cout << myqp << std::endl;

	myqp.set_x(x_global);
	myqp.set_x0(x0_global);
	myqp.set_b(b_global);
	myqp.set_A(A_global);
	
	std::cout << myqp << std::endl;

//pascinference::DEBUG_MODE = 0;	



	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

