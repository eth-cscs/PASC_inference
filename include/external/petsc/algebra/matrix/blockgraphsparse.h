#ifndef PASC_BLOCKGRAPHSPARSEMATRIX_PETSC_H
#define	PASC_BLOCKGRAPHSPARSEMATRIX_PETSC_H

#include "algebra/vector/generalvector.h"
#include "algebra/matrix/blockgraphsparse.h"

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
