#ifndef PASC_PETSCVECTOR_BLOCKGRAPHSPARSEMATRIX_H
#define	PASC_PETSCVECTOR_BLOCKGRAPHSPARSEMATRIX_H

#include "general/algebra/matrix/blockgraphsparse.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/matrix/generalmatrix.h"

namespace pascinference {
namespace algebra {

template<> BlockGraphSparseMatrix<PetscVector>::BlockGraphSparseMatrix(Decomposition<PetscVector> &new_decomposition, double alpha, GeneralVector<PetscVector> *new_coeffs);
template<> BlockGraphSparseMatrix<PetscVector>::~BlockGraphSparseMatrix();
template<> void BlockGraphSparseMatrix<PetscVector>::printcontent(ConsoleOutput &output) const;
template<> void BlockGraphSparseMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const;
template<> Mat BlockGraphSparseMatrix<PetscVector>::get_petscmatrix() const;

}
} /* end of namespace */

#endif
