
#include "algebra/vector/generalvector.h"

namespace pascinference {
namespace algebra {

template<>
void GeneralVector<PetscVector>::set_random() { 
	PetscRandom rnd;
	
	/* prepare random generator */
//	TRYCXX( PetscRandomCreate(PETSC_COMM_WORLD,&(this->rnd)) );
	TRYCXX( PetscRandomCreate(PETSC_COMM_SELF,&rnd) );

	TRYCXX( PetscRandomSetType(rnd,PETSCRAND) );
	TRYCXX( PetscRandomSetFromOptions(rnd) );

	TRYCXX( PetscRandomSetSeed(rnd,13) );

	Vec vec = this->get_vector();

	/* generate random data to gamma */
	TRYCXX( VecSetRandom(vec, rnd) );

	/* destroy the random generator */
	TRYCXX( PetscRandomDestroy(&rnd) );

	this->valuesUpdate();
}

template<>
void GeneralVector<PetscVector>::set_random2() { 
	PetscRandom rnd;

	/* prepare random generator */
	TRYCXX( PetscRandomCreate(PETSC_COMM_SELF,&rnd) );

	TRYCXX( PetscRandomSetType(rnd,PETSCRAND) );
	TRYCXX( PetscRandomSetFromOptions(rnd) );

	TRYCXX( PetscRandomSetSeed(rnd,13) );

	Vec vec = this->get_vector();

	/* random generator based on one-processor generator */
	Vec random_vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->size(), &random_vec) );
	TRYCXX( VecSetRandom(random_vec, rnd) );
	TRYCXX( VecAssemblyBegin(random_vec) );
	TRYCXX( VecAssemblyEnd(random_vec) );

	int local_begin;
	TRYCXX( VecGetOwnershipRange(vec, &local_begin, NULL) );
	
	IS local_is;
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, this->local_size(), local_begin, 1, &local_is) );

	Vec local_vec, local_random_vec;
	TRYCXX( VecGetSubVector(vec, local_is, &local_vec) );
	TRYCXX( VecGetSubVector(random_vec, local_is, &local_random_vec) );

	TRYCXX( VecCopy(local_random_vec, local_vec) );

	TRYCXX( VecRestoreSubVector(vec, local_is, &local_vec) );
	TRYCXX( VecRestoreSubVector(random_vec, local_is, &local_random_vec) );

	TRYCXX( ISDestroy(&local_is) );
	TRYCXX( VecDestroy(&random_vec) );

	/* destroy the random generator */
	TRYCXX( PetscRandomDestroy(&rnd) );

	this->valuesUpdate();
}


template<>
std::string GeneralVector<PetscVector>::get_name() {
	return "PetscVector";
}


}
} /* end of namespace */
