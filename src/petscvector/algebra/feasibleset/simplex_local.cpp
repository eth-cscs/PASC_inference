#include "external/petscvector/algebra/feasibleset/simplex_local.h"

namespace pascinference {
namespace algebra {

template<>
SimplexFeasibleSet_Local<PetscVector>::SimplexFeasibleSet_Local(int T, int K){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->T = T;
	this->K = K;

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();
	
	#ifdef USE_CUDA
		externalcontent->cuda_create(T,K);
	#endif

	LOG_FUNC_END
}

/* general destructor */
template<>
SimplexFeasibleSet_Local<PetscVector>::~SimplexFeasibleSet_Local(){
	LOG_FUNC_BEGIN
	
	#ifdef USE_CUDA
		externalcontent->cuda_destroy();
	#endif	
	
	LOG_FUNC_END	
}


template<>
void SimplexFeasibleSet_Local<PetscVector>::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	/* get local array */
	double *x_arr;
	
	#ifdef USE_CUDA
		TRYCXX( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );

		/* use kernel to compute projection */
		//TODO: here should be actually the comparison of Vec type! not simple use_gpu
		externalcontent->cuda_project(x_arr,T,K);

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

