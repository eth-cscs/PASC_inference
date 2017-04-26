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



}
} /* end of namespace */

