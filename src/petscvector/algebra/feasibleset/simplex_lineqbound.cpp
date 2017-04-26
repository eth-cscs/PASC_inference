#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"

namespace pascinference {
namespace algebra {

template<>
SimplexFeasibleSet_LinEqBound<PetscVector>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K){
	LOG_FUNC_BEGIN

	this->T = T;
	this->Tlocal = Tlocal;
	this->K = K;

	/* create global matrix from local blocks */
	TRYCXX( MatCreate(PETSC_COMM_WORLD, &B) );
	TRYCXX( MatSetSizes(B,this->Tlocal,this->K*this->Tlocal,this->T,this->K*this->T) );
	TRYCXX( MatSetFromOptions(B) );
	TRYCXX( MatMPIAIJSetPreallocation(B,this->K,NULL,this->K,NULL) );
	TRYCXX( MatSeqAIJSetPreallocation(B,this->K,NULL) );

	TRYCXX( MatAssemblyBegin(B,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(B,MAT_FLUSH_ASSEMBLY) );

	int row_begin;
	int col_begin;
	TRYCXX( MatGetOwnershipRange(B, &row_begin, NULL) );
	TRYCXX( MatGetOwnershipRangeColumn(B, &col_begin, NULL) );

	/* set local content */
	double value = 1.0;///(double)(this->K*this->T);
	for(int t=0; t<this->Tlocal;t++){
		for(int k=0;k<this->K;k++){
			TRYCXX( MatSetValue(B, row_begin + t, col_begin + t*this->K+k, value, INSERT_VALUES) );
		}
	}

	TRYCXX( MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)B,"equality constraint") );

	/* create RHS vector of equality constraints */
	TRYCXX( MatCreateVecs(B, &lb, &c) );
	TRYCXX( VecSet(c,value) );
	TRYCXX( VecAssemblyBegin(c) );
	TRYCXX( VecAssemblyEnd(c) );
	TRYCXX( PetscObjectSetName((PetscObject)c,"RHS eq") );

	/* create RHS vector of lower bounds */
	TRYCXX( VecSet(lb,0.0) );
	TRYCXX( VecAssemblyBegin(lb) );
	TRYCXX( VecAssemblyEnd(lb) );
	TRYCXX( PetscObjectSetName((PetscObject)lb,"lower bounds") );

	LOG_FUNC_END
}

template<>
Mat SimplexFeasibleSet_LinEqBound<PetscVector>::get_B() const {
	return this->B;
}

template<>
Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_c() const {
	return this->c;
}

template<>
Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_lb() const {
	return this->lb;
}



}
} /* end of namespace */

