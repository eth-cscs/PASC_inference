#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"

namespace pascinference {
namespace algebra {

template<>
SimplexFeasibleSet_LinEqBound<PetscVector>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K){
	LOG_FUNC_BEGIN

	this->T = T;
	this->Tlocal = Tlocal;
	this->K = K;

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();

	/* create global matrix from local blocks */
	TRYCXX( MatCreate(PETSC_COMM_WORLD, &(externalcontent->B)) );
	TRYCXX( MatSetSizes(externalcontent->B,this->Tlocal,this->K*this->Tlocal,this->T,this->K*this->T) );
	TRYCXX( MatSetFromOptions(externalcontent->B) );
	TRYCXX( MatMPIAIJSetPreallocation(externalcontent->B,this->K,NULL,this->K,NULL) );
	TRYCXX( MatSeqAIJSetPreallocation(externalcontent->B,this->K,NULL) );

	TRYCXX( MatAssemblyBegin(externalcontent->B,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(externalcontent->B,MAT_FLUSH_ASSEMBLY) );

	int row_begin;
	int col_begin;
	TRYCXX( MatGetOwnershipRange(externalcontent->B, &row_begin, NULL) );
	TRYCXX( MatGetOwnershipRangeColumn(externalcontent->B, &col_begin, NULL) );

	/* set local content */
	double value = 1.0;///(double)(this->K*this->T);
	for(int t=0; t<this->Tlocal;t++){
		for(int k=0;k<this->K;k++){
			TRYCXX( MatSetValue(externalcontent->B, row_begin + t, col_begin + t*this->K+k, value, INSERT_VALUES) );
		}
	}

	TRYCXX( MatAssemblyBegin(externalcontent->B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(externalcontent->B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)(externalcontent->B),"equality constraint") );

	/* create RHS vector of equality constraints */
	TRYCXX( MatCreateVecs(externalcontent->B, &(externalcontent->lb), &(externalcontent->c)) );
	TRYCXX( VecSet(externalcontent->c,value) );
	TRYCXX( VecAssemblyBegin(externalcontent->c) );
	TRYCXX( VecAssemblyEnd(externalcontent->c) );
	TRYCXX( PetscObjectSetName((PetscObject)(externalcontent->c),"RHS eq") );

	/* create RHS vector of lower bounds */
	TRYCXX( VecSet(externalcontent->lb,0.0) );
	TRYCXX( VecAssemblyBegin(externalcontent->lb) );
	TRYCXX( VecAssemblyEnd(externalcontent->lb) );
	TRYCXX( PetscObjectSetName((PetscObject)(externalcontent->lb),"lower bounds") );

	LOG_FUNC_END
}

template<> void SimplexFeasibleSet_LinEqBound<PetscVector>::project(GeneralVector<PetscVector> &x){
	LOG_FUNC_BEGIN

	Vec x_Vec = x.get_vector();
	
	/* get local size */
	int local_size;
	TRYCXX( VecGetLocalSize(x_Vec, &local_size) );
	
	/* projection onto box */
	double *x_arr;
	TRYCXX( VecGetArray(x_Vec, &x_arr) );
	for(int i=0;i<local_size;i++){
		if(x_arr[i] > 1.0){
			x_arr[i] = 1.0;
		}
		if(x_arr[i] < 0.0){
			x_arr[i] = 0.0;
		}
	}
	TRYCXX( VecRestoreArray(x_Vec, &x_arr) );
	
	LOG_FUNC_END
}

template<> SimplexFeasibleSet_LinEqBound<PetscVector>::ExternalContent * SimplexFeasibleSet_LinEqBound<PetscVector>::get_externalcontent() const {
	return this->externalcontent;
}


}
} /* end of namespace */

