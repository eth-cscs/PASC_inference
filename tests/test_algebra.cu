#include "pascinference.h"


#include "matrix/laplacefull.h"


using namespace pascinference;


int main( int argc, char *argv[] )
{
		
	Initialize(argc, argv); // TODO: load parameters from console input

	/* say hello */	
	Message("- start program");

	int N = 10;
	GlobalVector vg(N);
	HostVector vh(N);

	std::cout << "v_global: " << vg << std::endl;
	std::cout << "v_host: " << vh << std::endl;
	
	LaplaceFullMatrix<GlobalVector> Ag(vg);
	LaplaceFullMatrix<HostVector> Ah(vh);

	std::cout << "A_global: " << Ag << std::endl;
	std::cout << "A_host: " << Ah << std::endl;


	/* say bye */	
	Message("- end program");
	Finalize();
	return 0;
}

