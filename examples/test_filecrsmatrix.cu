#include "pascinference.h"
#include "matrix/filecrs.h"

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

	int N = 3;

	#ifdef USE_PETSCVECTOR 
		Initialize(argc, argv); // TODO: load parameters from console input
		petscvector::PETSC_INITIALIZED = true;

		std::cout << "-------------------- TEST GLOBAL --------------------" << std::endl;

		GeneralVector<Global> vg(N);
		vg(gall) = 3.3;
		std::cout << "v_global: " << vg << std::endl;

		FileCRSMatrix<Global> Ag(vg,"../three.bin");

		std::cout << "A_global: " << Ag << std::endl;
		GeneralVector<Global> Avg(N);

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

		FileCRSMatrix<Host> Ah(vh,"../three.bin");

		std::cout << "A_host:  " << Ah << std::endl;
		GeneralVector<Host> Avh(N);

		Avh = Ah*vh; 
		std::cout << "Av_host: " << Avh << std::endl;


		std::cout << "-------------------- TEST DEVICe --------------------" << std::endl;

		GeneralVector<Device> vd(N);
		vd(gall) = 3.3;
		std::cout << "v_device:  " << vd << std::endl;

		FileCRSMatrix<Device> Ad(vd,"../three.bin");

		std::cout << "A_device:  " << Ad << std::endl;
		GeneralVector<Device> Avd(N);

		Avd = Ad*vd; 
		std::cout << "Av_device: " << Avd << std::endl;

	#endif

	/* say bye */	
	Message("- end program");
	return 0;
}

