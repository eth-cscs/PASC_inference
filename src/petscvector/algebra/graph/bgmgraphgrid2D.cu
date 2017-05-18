#include "external/petscvector/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace algebra {

void BGMGraphGrid2D<PetscVector>::ExternalContent::process_grid_cuda(){
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
} /* end of namespace */
