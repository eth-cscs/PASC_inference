#include "external/petscvector/data/entropydata.h"

namespace pascinference {
namespace data {

/* constructor */
template<>
EntropyData<PetscVector>::EntropyData(int T, int xdim, int K, int Km){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->lambda = NULL;
	this->x = NULL;
	this->gamma = NULL;

	this->T = T;
	this->xdim = xdim;
	this->K = K;
	this->Km = Km;

	this->number_of_moments = compute_number_of_moments(this->xdim, this->Km);

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();

	/* prepare matrix D */
	this->prepare_matrix_D();

	LOG_FUNC_END
}

template<>
void EntropyData<PetscVector>::prepare_matrix_D() {
	LOG_FUNC_BEGIN

	Vec matrix_D_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&matrix_D_Vec) );
	TRYCXX( VecSetType(matrix_D_Vec, VECSEQ));
	TRYCXX( VecSetSizes(matrix_D_Vec, get_xdim()*get_number_of_moments(), PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(matrix_D_Vec) );

	double *matrix_D_arr;
	TRYCXX( VecGetArray(matrix_D_Vec, &matrix_D_arr) );

	/* run outer for cycle and call recursion */
	int idx = 0;
	externalcontent->prepare_matrix_D_recursion(matrix_D_arr, &idx, get_Km(), 0, get_Km(), get_xdim());

	TRYCXX( VecRestoreArray(matrix_D_Vec, &matrix_D_arr) );

	this->matrix_D = new GeneralVector<PetscVector>(matrix_D_Vec);

	LOG_FUNC_END
}

template<>
void EntropyData<PetscVector>::print_matrix_D(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	double *matrix_D_arr;
	TRYCXX( VecGetArray(this->matrix_D->get_vector(), &matrix_D_arr) );

	for(int i_moment=0;i_moment < get_number_of_moments(); i_moment++){
		for(int i_xdim=0;i_xdim < get_xdim(); i_xdim++){
			output << matrix_D_arr[ i_moment*get_xdim() + i_xdim ];
			if(i_xdim < get_xdim()-1) output << ", ";
		}
		output << std::endl;
	}

	TRYCXX( VecRestoreArray(this->matrix_D->get_vector(), &matrix_D_arr) );

	LOG_FUNC_END
}



template<> EntropyData<PetscVector>::ExternalContent * EntropyData<PetscVector>::get_externalcontent() const {
	return externalcontent;
}

void EntropyData<PetscVector>::ExternalContent::prepare_matrix_D_recursion(double *values, int *idx, int top_i, int level, int Km, int xdim){

	for(int i=0; i<=top_i; i++){ /* FOR cycle on this level */
		/* copy values from upper row to new row
			(represents constant iterators of outer for cycles) */
		if(i>0){
			for(int previous_level=0;previous_level <= level-1; previous_level++){
				double temp = get_D_value(values, *idx - 1, previous_level, xdim);
				set_D_value(values, temp, *idx, previous_level, xdim);
			}
		}

		/* add my own value */
		set_D_value(values, i, *idx, level, xdim);

		if(level < xdim - 1){
			/* go deeper with recursion */
			prepare_matrix_D_recursion(values, idx, top_i - i, level+1, Km, xdim);
		} else {
			/* last level wrote the last value into this row */
			*idx += 1;
		}
	}

}

void EntropyData<PetscVector>::ExternalContent::set_D_value(double *values, double value, int row, int col, int ncols){
	coutMaster << "value: " << value << std::endl;
	coutMaster << "row:   " << row << std::endl;
	coutMaster << "col:   " << col << std::endl;
	coutMaster << "ncols: " << ncols << std::endl;
	coutMaster << "index: " << row*ncols + col << std::endl;

	values[row*ncols + col] = value;
}

double EntropyData<PetscVector>::ExternalContent::get_D_value(double *values, int row, int col, int ncols) const {
	return values[row*ncols + col];
}


}
} /* end namespace */

