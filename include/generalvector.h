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

	/* deal with all */
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

	} gall; /* sorry for gall, but minlin is using global all and I don't know how to retype it. */ // TODO: deal with all


	/* general vector class - take original class and add multiplication with GeneralMatrix */
	template<class VectorBase>
	class GeneralVector : public VectorBase {

		public:
			/* constructors */
			GeneralVector(): VectorBase() {}
			template<class ArgType> GeneralVector(ArgType arg): VectorBase(arg) {}
			
			/* matrix-vector multiplication with General matrix */
			GeneralVector<VectorBase> &operator=(GeneralMatrixRHS<VectorBase> rhs){
				rhs.matmult(*this);
				return *this;
			}

			/* fill with random values */
			void set_random();

	};

}

// TODO: move to impl
namespace pascinference {

/* ------- PETSCVECTOR ------- */	
#ifdef USE_PETSCVECTOR
typedef petscvector::PetscVector PetscVector;

template<>
void GeneralVector<PetscVector>::set_random() { 
	PetscRandom rnd;
	Vec vec = this->get_vector();
	
	/* prepare random generator */
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );

	TRY( PetscRandomSetSeed(rnd,13) );

	/* generate random data to gamma */
	TRY( VecSetRandom(vec, rnd) );

	/* destroy the random generator */
	TRY( PetscRandomDestroy(&rnd) );

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
