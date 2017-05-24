#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {
	
void BGMGraph<PetscVector>::ExternalContent::cuda_destroy(){
	LOG_FUNC_BEGIN

	gpuErrchk( cudaFree(neighbor_nmbs_gpu) );
	for(int i=0;i<n;i++){
		gpuErrchk( cudaFree(neighbor_ids_cpugpu[i]) );
	}
	free(neighbor_ids_cpugpu);
	gpuErrchk( cudaFree(neighbor_ids_gpu) );

	LOG_FUNC_END
}

void BGMGraph<PetscVector>::ExternalContent::cuda_process(int *neighbor_nmbs, int **neighbor_ids){
	LOG_FUNC_BEGIN

	/* copy data to gpu */
	gpuErrchk( cudaMalloc((void **)&neighbor_nmbs_gpu, n*sizeof(int)) );
	gpuErrchk( cudaMemcpy( neighbor_nmbs_gpu, neighbor_nmbs, n*sizeof(int), cudaMemcpyHostToDevice) );
		
	/* allocate pointers on CPU */
	neighbor_ids_cpugpu = (int**)malloc(n*sizeof(int*));
		
	for(int i=0;i<n;i++){
		int mysize = neighbor_nmbs[i];
		
		gpuErrchk( cudaMalloc((void **)&(neighbor_ids_cpugpu[i]), mysize*sizeof(int)) );
		gpuErrchk( cudaMemcpy( neighbor_ids_cpugpu[i], neighbor_ids[i], mysize*sizeof(int), cudaMemcpyHostToDevice) );
	}

	/* copy pointers to arrays from CPU to GPU */
	gpuErrchk( cudaMalloc((void **)&neighbor_ids_gpu, n*sizeof(int*)) );
	gpuErrchk( cudaMemcpy( neighbor_ids_gpu, neighbor_ids_cpugpu, n*sizeof(int*), cudaMemcpyHostToDevice) );

	gpuErrchk( cudaDeviceSynchronize() );

	LOG_FUNC_END
}

}
}
