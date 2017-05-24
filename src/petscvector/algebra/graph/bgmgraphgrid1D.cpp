#include "external/petscvector/algebra/graph/bgmgraphgrid1D.h"

namespace pascinference {
namespace algebra {

template<>
BGMGraphGrid1D<PetscVector>::BGMGraphGrid1D(int width) : BGMGraph<PetscVector>(){
	LOG_FUNC_BEGIN

	this->width = width;

	this->dim = 2;
	this->n = width;
	
	/* fill coordinates */
	Vec coordinates_Vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );
	
	double *coordinates_arr;
	TRYCXX( VecGetArray(coordinates_Vec, &coordinates_arr) );
	for(int i=0;i<width;i++){
		coordinates_arr[i] = i;
		coordinates_arr[i + this->n] = 0;
	}
	TRYCXX( VecRestoreArray(coordinates_Vec, &coordinates_arr) );
	
	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	this->processed = false;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid1D<PetscVector>::process_grid(){
	this->threshold = 1.1;
	this->m = width-1;
	this->m_max = 2;

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

	for(int i=0;i<width;i++){
		int idx = i;

		/* compute number of neighbors */
		int nmb = 0;
		if(i>0){
			nmb+=1;				
		}
		if(i<width-1){
			nmb+=1;				
		}
		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));
			
		/* fill neighbors */
		nmb = 0;
		if(i>0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb++;
		}
		if(i<width-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb++;
		}
	}

	#ifdef USE_CUDA
		externalcontent->process_grid_cuda();
	#endif
	
	this->processed = true;
}

template<> BGMGraphGrid1D<PetscVector>::ExternalContent * BGMGraphGrid1D<PetscVector>::get_externalcontent() const {
	return externalcontent;
}


}
} /* end of namespace */

