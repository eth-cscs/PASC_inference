#include "pascinference.h"
#include "matrix/blockdiaglaplace_explicit.h"

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

	int T = 5; /* size of one block */
	int K = 3; /* number of block */
	int N = T*K;

	std::cout << "T = " << T << std::endl;
	std::cout << "K = " << K << std::endl;
	std::cout << "N = " << N << std::endl;

	#ifdef USE_PETSCVECTOR
		Initialize(argc, argv); // TODO: load parameters from console input
		petscvector::PETSC_INITIALIZED = true;

		std::cout << "-------------------- TEST GLOBAL --------------------" << std::endl;

		GeneralVector<Global> vg(N);
		vg(gall) = 3.3;
		std::cout << "v_global: " << vg << std::endl;

		BlockDiagLaplaceExplicitMatrix<Global> Ag(vg,K, 5.5);
		
		std::cout << "A_global: " << Ag << std::endl;
		GeneralVector<Global> Avg(vg);

		Avg = Ag*vg; 
		std::cout << "Av_global: " << Avg << std::endl;

		petscvector::PETSC_INITIALIZED = false;
		Finalize();

	#endif


	#ifdef USE_MINLIN
		std::cout << "-------------------- TEST HOST  --------------------" << std::endl;

		GeneralVector<Host> vh(N);
		vh(gall) = 3.3;
		std::cout << "v_host:  " << vh << std::endl;

		BlockDiagLaplaceExplicitMatrix<Host> Ah(vh,K);

		std::cout << "A_host:  " << Ah << std::endl;
		GeneralVector<Host> Avh(vh);

		Avh = Ah*vh; 
		std::cout << "Av_host: " << Avh << std::endl;



		std::cout << "-------------------- TEST DEVICe --------------------" << std::endl;

		GeneralVector<Device> vd(N);
		vd(gall) = 3.3;
		std::cout << "v_device:  " << vd << std::endl;

		BlockDiagLaplaceExplicitMatrix<Device> Ad(vd,K);

		std::cout << "A_device:  " << Ad << std::endl;
		GeneralVector<Device> Avd(vd);

		Avd = Ad*vd; 
		std::cout << "Av_device: " << Avd << std::endl;

	#endif

	/* say bye */	
	Message("- end program");
	return 0;
}

