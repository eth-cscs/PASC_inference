#include "pascinference.h"
#include "matrix/laplace_explicit.h"

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


int main( int argc, char *argv[] )
{
		
	/* say hello */	
	Message("- start program");

	int N = 5;

	#ifdef USE_PETSCVECTOR
		Message("------------- PETSC TEST --------------");
		Initialize(argc, argv);
		petscvector::PETSC_INITIALIZED = true;

		GeneralVector<Global> vg(N);
		vg(gall) = 3.3;
		std::cout << "v_global: " << vg << std::endl;
	
		LaplaceExplicitMatrix<Global> Ag(vg);
		std::cout << "A_global: " << Ag << std::endl;

		GeneralVector<Global> Avg(N);
		Avg = Ag*vg; 
		std::cout << "Av_global: " << Avg << std::endl;

		petscvector::PETSC_INITIALIZED = false;
		Finalize();

	#endif

	#ifdef USE_MINLIN
		Message("------------- MINLIN Host TEST --------------");
		GeneralVector<Host> vh(N);
		vh(gall) = 3.3;
		std::cout << "v_host: " << vh << std::endl;
	
		LaplaceExplicitMatrix<Host> Ah(vh);
		std::cout << "A_host: " << Ah << std::endl;

		GeneralVector<Host> Avh(N);
		Avh = Ah*vh; 
		std::cout << "Av_host: " << Avh << std::endl;
	#endif

	#ifdef USE_MINLIN
		Message("------------- MINLIN Device TEST --------------");
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
	

	return 0;
}

