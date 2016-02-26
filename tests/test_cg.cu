#include "pascinference.h"
#include "solver/cg.h"
#include "matrix/laplace_explicit.h"


using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
typedef minlin::threx::HostVector<double> Host;
typedef minlin::threx::DeviceVector<double> Device;

extern bool petscvector::PETSC_INITIALIZED;


int main( int argc, char *argv[] )
{
		
	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;
	
	/* say hello */	
	Message("- start program");

	/* dimension of the problem */
	int N = 1000;

/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare solution vector */
	Vector<Global> x_global(N);
	x_global(gall) = 0.0;

	/* prepare initial approximation */
	Vector<Global> x0_global(x_global);
	x0_global(gall) = 0.0;

	/* prepare right hand-side vector */
	Vector<Global> b_global(x_global);
	b_global(gall) = 1.0/N;

	/* prepare matrix */	
	LaplaceExplicitMatrix<Global> A_global(x_global);

	/* solve problem */
	x_global = cg(A_global, b_global, x0_global);
	
	
	/* print solution */
	std::cout << "x_global: " << x_global << std::endl;




	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

