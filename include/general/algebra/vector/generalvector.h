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

#include "general/algebra/matrix/generalmatrixrhs.h"
#include "general/algebra/arrayoperation.h"

//#include <iostream>
#include <string>

namespace pascinference {
namespace algebra {

	/** @class General_all_type
	*  @brief deal with vector(all)
	* 
	*  For each vector type define retype operator of our "all".
	* 
	*  @todo improve "gall", who to change it to "all"
	*/ 
	class General_all_type {
		public:
			#ifdef USE_PETSC
				/* convert to PetscVector all */
//TODO:
//				operator petscvector::petscvector_all_type() const { 
//					return petscvector::all;
//				}
			#endif

			#ifdef USE_MINLIN
				typedef minlin::detail::all_type *minlin_all_type;
				/* convert to minlin all */ 
				operator minlin_all_type() const {
					return minlin::all;
				}
			#endif

	};
	static General_all_type gall;


	/** @class GeneralVector
	*  @brief general vector class
	* 
	*  Extend the operations with original vector type.
	*  Add operator for multiplication with GeneralMatrix.
	* 
	*/ 
	template<class VectorBase>
	class GeneralVector : public VectorBase {
		public:
			/** @brief call original constructor without arguments
			*
			*/
			GeneralVector() : VectorBase() {
			}

			/** @brief call original constructor with one argument
			*
			* @todo for general number of arguments
			*/
			template<class ArgType> GeneralVector(ArgType arg) : VectorBase(arg) {
			}

			/** @brief call original constructor with two arguments
			*
			* @todo for general number of arguments
			*/
			template<class ArgType1,class ArgType2> GeneralVector(ArgType1 arg1, ArgType2 arg2) : VectorBase(arg1,arg2){
			}
			
			/** @brief matrix vector multiplication
			*
			*  Set the values of vector equal to the values obtained from right hand-side 
			*  vector A*v.
			*  \f[ \mathrm{this}~ = \underbrace{Av}_{= \mathrm{rhs}} \f]
			*  
			*  @param rhs right hand-size vector from A*v
			*/
			GeneralVector<VectorBase> &operator=(GeneralMatrixRHS<VectorBase> rhs);

			/** @brief set random values
			*
			*/
			virtual void set_random();

			/** @brief set random values
			*
			*/
			virtual void set_random2();
			
			/** @brief get the name of the type
			*
			*/
			static std::string get_name();
			
	};

}
}

/* ------ IMPLEMENTATION ------ */
namespace pascinference {
namespace algebra {

template<class VectorBase>
GeneralVector<VectorBase> & GeneralVector<VectorBase>::operator=(GeneralMatrixRHS<VectorBase> rhs){
	rhs.matmult(*this);
	return *this;
}

template<class VectorBase>
void GeneralVector<VectorBase>::set_random(){
}

template<class VectorBase>
void GeneralVector<VectorBase>::set_random2(){
}
			
template<class VectorBase>
std::string GeneralVector<VectorBase>::get_name(){
	return "GeneralVector";
}

			
}
}



// TODO: move to impl
namespace pascinference {
namespace algebra {

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

}
} /* end of namespace */

#endif
