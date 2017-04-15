#ifndef PASC_GENERALVECTOR_PETSC_H
#define	PASC_GENERALVECTOR_PETSC_H

#include "algebra/vector/generalvector.h"
#include "external/petsc/algebra/vector/petscvector.h"

typedef petscvector::PetscVector PetscVector;

namespace pascinference {

namespace algebra {

template class GeneralVector<PetscVector>;
template class GeneralMatrixRHS<PetscVector>;

//template<> GeneralVector<PetscVector>::GeneralVector();
//extern template GeneralVector<PetscVector>::GeneralVector(int n);
//extern template GeneralVector<PetscVector>::GeneralVector(double *values, int n);
//template<> GeneralVector<PetscVector>::GeneralVector(const PetscVector &vec1);
//template<> template<> GeneralVector<PetscVector>::GeneralVector<const Vec &>(const Vec &new_inner_vector);
//extern template GeneralVector<PetscVector>::GeneralVector(const PetscVectorWrapperComb &comb);


//static GeneralMatrixRHS<PetscVector> bbb;

//template<> class GeneralVector<PetscVector>;
//extern GeneralVector<PetscVector>::GeneralVector(Vec myvec);

//template<> void GeneralVector<PetscVector>::set_random();
//template<> void GeneralVector<PetscVector>::set_random2();
//template<> std::string GeneralVector<PetscVector>::get_name();


//extern template class GeneralMatrixRHS<PetscVector>;
//extern template class GeneralVector<PetscVector>;

//template<class VectorBase> class GeneralVector;
//extern template<> void GeneralVector<PetscVector>::set_random();
//extern template void GeneralVector<PetscVector>::set_random2();
//extern template std::string GeneralVector<PetscVector>::get_name();

//template<>
//extern void GeneralVector<PetscVector>::set_random();

}
} /* end of namespace */

#endif
