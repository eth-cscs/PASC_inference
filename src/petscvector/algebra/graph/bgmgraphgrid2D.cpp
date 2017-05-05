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

//	#pragma omp parallel for
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

}
} /* end of namespace */
