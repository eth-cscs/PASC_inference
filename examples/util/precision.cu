/** @file precision.cu
 *  @brief test the precision of stored values (float/double?)
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <math.h>

using namespace pascinference;

#ifdef USE_CUDA
__global__ void print_precision(){
	double a;
	for(int i=0; i < 25;i++){
		a = 1.0 + pow(10.0,-i);
		printf(" 1+10^{- %d} : %.30f\n", i, a );
	}
}

#endif


int main( int argc, char *argv[] )
{
	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 
	
	coutMaster << std::setprecision(30);
	coutMaster << "CPU:" << std::endl;

	double a;
	for(int i=0; i < 25;i++){
		a = 1.0 + pow(10.0,-i);
		coutMaster << " 1+10^{-" << std::setw(3) << i << "}: " << a << std::endl;
	}


#ifdef USE_CUDA
	coutMaster << std::endl;
	coutMaster << "GPU:" << std::endl;

	print_precision<<<1, 1>>>();
	gpuErrchk( cudaDeviceSynchronize() );

#endif

	Finalize();
	
	return 0;
}

