/** @file test_algebra.cpp
 *  @brief test matrix-vector multiplication
 *
 *  Create laplace explicit matrix and vector. Multiply these objects on different architectures.
 *
 *  @author Lukas Pospisil
 */


#include "pascinference.h"
#include "matrix/laplace_explicit.h"

using namespace pascinference;

/* set what is what ( which type of vector to use where) */
#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector Global;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> Host;
	typedef minlin::threx::DeviceVector<double> Device;
#endif


int main( int argc, char *argv[] )
{
		
	/* say hello */	
	Message("- start program");
	Initialize(argc, argv);

	int N = 5;

	/* ------------- PETSC TEST -------------- */
	#ifdef USE_PETSCVECTOR

		GeneralVector<Global> vg(N);
		vg(gall) = 3.3;
		std::cout << "v_global: " << vg << std::endl;
	
		LaplaceExplicitMatrix<Global> Ag(vg);
		std::cout << "A_global: " << Ag << std::endl;

		GeneralVector<Global> Avg(N);
		Avg = Ag*vg; 
		std::cout << "Av_global: " << Avg << std::endl;

		petscvector::PETSC_INITIALIZED = false;

	#endif

	/* ------------- MINLIN TEST -------------- */
	#ifdef USE_MINLIN
		GeneralVector<Host> vh(N);
		vh(gall) = 3.3;
		std::cout << "v_host: " << vh << std::endl;
	
		LaplaceExplicitMatrix<Host> Ah(vh);
		std::cout << "A_host: " << Ah << std::endl;

		GeneralVector<Host> Avh(N);
		Avh = Ah*vh; 
		std::cout << "Av_host: " << Avh << std::endl;


		GeneralVector<Device> vd(N);
		vd(gall) = 3.3;
		std::cout << "v_device: " << vd << std::endl;
	
		LaplaceExplicitMatrix<Device> Ad(vd);
		std::cout << "A_device: " << Ad << std::endl;

		GeneralVector<Device> Avd(N);
		Avd = Ad*vd; 
		std::cout << "Av_device: " << Avd << std::endl;

	#endif


	/* say bye */	
	Message("- end program");
	Finalize();

	return 0;
}

