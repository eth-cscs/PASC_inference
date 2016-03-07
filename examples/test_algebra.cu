#include "pascinference.h"


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

	int N = 5;
	GeneralVector<Global> vg(N);
	GeneralVector<Host> vh(N);

	vg(0) = 1.0;
	vg(1) = 2.0;
	vg(2) = 3.0;
	vg(3) = 4.0;
	vg(4) = 5.0;
	
	vg(gall) = 3.3;
	vh(gall) = 2.3;

	std::cout << "v_global: " << vg << std::endl;
	std::cout << "v_host:   " << vh << std::endl;
	
	LaplaceExplicitMatrix<Global> Ag(vg);
	LaplaceExplicitMatrix<Host> Ah(vh);

	std::cout << "A_global: " << Ag << std::endl;
	std::cout << "A_host: " << Ah << std::endl;

	GeneralVector<Global> Avg(N);
	GeneralVector<Host> Avh(N);
	Avg = Ag*vg; 
	Avh = Ah*vh; 

	std::cout << "Av_global: " << Avg << std::endl;
	std::cout << "Av_host: " << Avh << std::endl;


	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

