#include "algebra/feasibleset/simplex_lineqbound.h"

#ifndef USE_PETSC
 #error 'SIMPLEX_LINEQBOUNDFEASIBLESET is for PETSC'
#endif

namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K){
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

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::~SimplexFeasibleSet_LinEqBound(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet_LinEqBound<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nmb of subsets:     " << T << std::endl;
	output <<  " - size of subset:     " << K << std::endl;


	LOG_FUNC_END
}

template<class VectorBase>
std::string SimplexFeasibleSet_LinEqBound<VectorBase>::get_name() const {
	return "SimplexFeasibleSet_LinEqBound";
}

template<>
void SimplexFeasibleSet_LinEqBound<PetscVector>::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	//TODO: give error, the projection is not defined for this type of feasible set
	TRYCXX( PetscBarrier(NULL) );

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
} /* end namespace */
