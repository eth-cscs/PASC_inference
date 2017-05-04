#include "external/petscvector/algebra/feasibleset/simplex_local.h"

namespace pascinference {
namespace algebra {

template<>
void SimplexFeasibleSet_Local<PetscVector>::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	/* get local array */
	double *x_arr;
	
	#ifdef USE_CUDA
		TRYCXX( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );

		/* use kernel to compute projection */
		//TODO: here should be actually the comparison of Vec type! not simple use_gpu
		kernel_project<<<gridSize, blockSize>>>(x_arr,x_sorted,T,K);
		gpuErrchk( cudaDeviceSynchronize() );
		MPI_Barrier( MPI_COMM_WORLD );

		TRYCXX( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	#else
		TRYCXX( VecGetArray(x.get_vector(),&x_arr) );
	
		/* use openmp */
//		#pragma omp parallel for
		for(int t=0;t<T;t++){
			project_sub(x_arr,t,T,K);
		}

		TRYCXX( VecRestoreArray(x.get_vector(),&x_arr) );
	#endif

	TRYCXX( PetscBarrier(NULL) );

	LOG_FUNC_END
}


}
} /* end of namespace */

