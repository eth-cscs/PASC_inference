/** @file generalvector.h
 *  @brief class for manipulation with vectors
 *
 *  Defines the parent class for manipulation with vectors.
 *  The class is inherited from other third-party implementations, i.e. MinLin or PetscVector.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_GENERALVECTOR_H
#define	PASC_GENERALVECTOR_H

namespace pascinference {
	/** @class General_all_type
	*  @brief deal with vector(all)
	* 
	*  For each vector type define retype operator of our "all".
	* 
	*  @todo improve "gall", who to change it to "all"
	*/ 
	class General_all_type {
		public:
			#ifdef USE_PETSCVECTOR
				/* convert to PetscVector all */
				operator petscvector::petscvector_all_type() const { 
					return petscvector::all;
				}
			#endif

			#ifdef USE_MINLIN
				typedef minlin::detail::all_type *minlin_all_type;
				/* convert to minlin all */ 
				operator minlin_all_type() const {
					return minlin::all;
				}
			#endif

	} gall;


	/** @class GeneralVector
	*  @brief general vector class
	* 
	*  Extend the operations with original vector type.
	*  Add operator for multiplication with GeneralMatrix.
	* 
	*/ 
	template<class VectorBase>
	class GeneralVector : public VectorBase {
		private:
			/* random vector generator */
			#ifdef USE_PETSCVECTOR
				PetscRandom rnd; /**< PETSc random generator */
			#endif

		public:
			/** @brief call original constructor without arguments
			*
			*/
			GeneralVector(): VectorBase() {
				#ifdef USE_PETSCVECTOR
					this->rnd=NULL;
				#endif
			}

			/** @brief call original constructor with one argument
			*
			* @todo for general number of arguments
			*/
			template<class ArgType> GeneralVector(ArgType arg): VectorBase(arg) {
				#ifdef USE_PETSCVECTOR
					this->rnd=NULL;
				#endif
			}

			/** @brief call original constructor with two arguments
			*
			* @todo for general number of arguments
			*/
			template<class ArgType1,class ArgType2> GeneralVector(ArgType1 arg1, ArgType2 arg2): VectorBase(arg1,arg2) {
				#ifdef USE_PETSCVECTOR
					this->rnd=NULL;
				#endif
			}
			
			/** @brief matrix vector multiplication
			*
			*  Set the values of vector equal to the values obtained from right hand-side 
			*  vector A*v.
			*  \f[ \mathrm{this}~ = \underbrace{Av}_{= \mathrm{rhs}} \f]
			*  
			*  @param rhs right hand-size vector from A*v
			*/
			GeneralVector<VectorBase> &operator=(GeneralMatrixRHS<VectorBase> rhs){
				rhs.matmult(*this);
				return *this;
			}

			/** @brief set random values
			*
			*/
			virtual void set_random();

	};

}

// TODO: move to impl
namespace pascinference {

/* ------- PETSCVECTOR ------- */	
#ifdef USE_PETSCVECTOR
typedef petscvector::PetscVector PetscVector;

template<>
void GeneralVector<PetscVector>::set_random() { 
	if(!this->rnd){
		/* prepare random generator */
//		TRY( PetscRandomCreate(PETSC_COMM_WORLD,&(this->rnd)) );
		TRY( PetscRandomCreate(PETSC_COMM_SELF,&(this->rnd)) );

		TRY( PetscRandomSetType(this->rnd,PETSCRAND) );
		TRY( PetscRandomSetFromOptions(this->rnd) );

		TRY( PetscRandomSetSeed(this->rnd,13) );
	}

	Vec vec = this->get_vector();

	/* generate random data to gamma */
//	TRY( VecSetRandom(vec, this->rnd) );

	/* random generator based on one-processor generator */
	Vec random_vec;
	TRY( VecCreateSeq(PETSC_COMM_SELF, this->size(), &random_vec) );
	TRY( VecSetRandom(random_vec, this->rnd) );
	TRY( VecAssemblyBegin(random_vec) );
	TRY( VecAssemblyEnd(random_vec) );

	int local_begin;
	TRY( VecGetOwnershipRange(vec, &local_begin, NULL) );
	
	IS local_is;
	TRY( ISCreateStride(PETSC_COMM_WORLD, this->local_size(), local_begin, 1, &local_is) );

	Vec local_vec, local_random_vec;
	TRY( VecGetSubVector(vec, local_is, &local_vec) );
	TRY( VecGetSubVector(random_vec, local_is, &local_random_vec) );

	TRY( VecCopy(local_random_vec, local_vec) );

	TRY( VecRestoreSubVector(vec, local_is, &local_vec) );
	TRY( VecRestoreSubVector(random_vec, local_is, &local_random_vec) );

	TRY( ISDestroy(&local_is) );
	TRY( VecDestroy(&random_vec) );

	/* destroy the random generator */
//	TRY( PetscRandomDestroy(&rnd) );

	this->valuesUpdate();
}

#endif

/* ------- MINLIN ------- */	
#ifdef USE_MINLIN
typedef minlin::threx::HostVector<double> MinlinHostVector;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;

template<>
void GeneralVector<MinlinHostVector>::set_random() { 
	int i;
	for(i=0;i<this->size();i++){
		(*this)(i) = std::rand()/(double)(RAND_MAX); /* generate random value */
	}	
}

template<>
void GeneralVector<MinlinDeviceVector>::set_random() { 
	int i;
	for(i=0;i<this->size();i++){ // TODO: fill on Host and then transfer to device
		(*this)(i) = std::rand()/(double)(RAND_MAX); /* generate random value */
	}	
}


#endif



} /* end of namespace petscvector */

#endif
