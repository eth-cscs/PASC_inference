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
	
	Vec x_Vec = x.get_vector();
	double *x_arr;
	
	#ifdef USE_CUDA
		/* use kernel to compute projection */
		TRYCXX( VecCUDAGetArrayReadWrite(x_Vec,&x_arr) );		

		externalcontent->cuda_project(x_Vec,T,K);

		TRYCXX( VecCUDARestoreArrayReadWrite(x_Vec,&x_arr) );		
	#else
		/* get local array */
		TRYCXX( VecGetArray(x_Vec,&x_arr) );
	
		for(int t=0;t<T;t++){
			project_sub(x_arr,t,T,K);
		}

		TRYCXX( VecRestoreArray(x_Vec,&x_arr) );
	#endif

	TRYCXX( PetscBarrier(NULL) );

	LOG_FUNC_END
}


}
} /* end of namespace */

