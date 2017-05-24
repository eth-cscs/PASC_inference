#include "external/petscvector/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace algebra {

template<>
BGMGraphGrid2D<PetscVector>::BGMGraphGrid2D(int width, int height) : BGMGraph<PetscVector>(){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;

	this->dim = 2;
	this->n = width*height;
	
	/* fill coordinates */
	Vec coordinates_Vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );
	
	double *coordinates_arr;
	TRYCXX( VecGetArray(coordinates_Vec, &coordinates_arr) );

	for(int idx=0;idx<width*height;idx++){
		int i = idx/(double)width; /* index of row */
		int j = idx - i*width; /* index of column */	

		coordinates_arr[idx] = j;
		coordinates_arr[idx + this->n] = i;
	}

	TRYCXX( VecRestoreArray(coordinates_Vec, &coordinates_arr) );
	
	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	this->processed = false;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid2D<PetscVector>::process_grid(){
	LOG_FUNC_BEGIN

	this->threshold = 1.1;
	this->m = height*(width-1) + width*(height-1);
	this->m_max = 4;

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int idx=0;idx<width*height;idx++){
		int i = idx/(double)width; /* index of row */
		int j = idx - i*width; /* index of column */	

		/* compute number of neighbors */
		int nmb = 0;
		if(j>0){
			nmb+=1;				
		}
		if(j<width-1){
			nmb+=1;				
		}
		if(i>0){
			nmb+=1;				
		}
		if(i<height-1){
			nmb+=1;				
		}
		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));
			
		/* fill neighbors */
		nmb = 0;
		if(j>0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb+=1;	
		}
		if(j<width-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb+=1;	
		}
		if(i>0){ /* down */
			this->neighbor_ids[idx][nmb] = idx-width;
			nmb+=1;	
		}
		if(i<height-1){ /* up */
			this->neighbor_ids[idx][nmb] = idx+width;
			nmb+=1;	
		}
	}

	#ifdef USE_CUDA
		externalcontent->process_grid_cuda(neighbor_nmbs, neighbor_ids);
	#endif
	
	this->processed = true;

	LOG_FUNC_END
}

template<> BGMGraphGrid2D<PetscVector>::ExternalContent * BGMGraphGrid2D<PetscVector>::get_externalcontent() const {
	return externalcontent;
}

}
} /* end of namespace */
