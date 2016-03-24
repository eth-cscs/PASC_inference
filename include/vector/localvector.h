#ifndef LOCALVECTOR_H
#define	LOCALVECTOR_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h" /* parent GeneralMatrix class */

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif

namespace pascinference {

/** @class LocalVector
*  @brief local vector class
* 
*  Provide operations with local vector. The vector could be initialized with given array of values.
* 
*/ 
template<class VectorBase>
class LocalVector: public GeneralVector<VectorBase> {
	private:
	
	public:
		/** @brief call original constructors
		*
		*/
		LocalVector(): GeneralVector<VectorBase>() {};

		/** @brief call original constructors
		*
		*/
		template<class ArgType> LocalVector(ArgType arg): GeneralVector<VectorBase>(arg) {};

		/** @brief call constructor with given array of values and size
		 * 
		 */
		LocalVector(double *values, int n){
				// TODO: give error, this function has to be specified for each type
		};

};

} /* end of namespace */


/* implementations */
//TODO: move to impls

/* -------------------------------- PETSC VECTOR -------------------------*/

namespace pascinference {

#ifdef USE_PETSCVECTOR

/* constructor from dimension creates local vector instead of MPI */
template<>
LocalVector<PetscVector>::LocalVector(double *values, int n) : GeneralVector<PetscVector>(values, n) {};

#endif

#ifdef USE_MINLIN

/* constructor from dimension creates local vector instead of MPI */
template<>
LocalVector<MinlinHostVector>::LocalVector(double *values, int n) : GeneralVector<MinlinHostVector>(n) {
	for(int i=0;i<n;i++){
		(*this)(i) = values[i];
	}
};

#endif


}

#endif
