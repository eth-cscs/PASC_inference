#ifndef PASC_PETSCVECTOR_BLOCKGRAPHSPARSETRMATRIX_H
#define	PASC_PETSCVECTOR_BLOCKGRAPHSPARSETRMATRIX_H

#include "general/algebra/matrix/blockgraphsparseTR.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/matrix/generalmatrix.h"
#include "external/petscvector/common/decomposition.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

template<> BlockGraphSparseTRMatrix<PetscVector>::BlockGraphSparseTRMatrix(Decomposition<PetscVector> &new_decomposition, double sigma, double alpha);
template<> BlockGraphSparseTRMatrix<PetscVector>::~BlockGraphSparseTRMatrix();
template<> void BlockGraphSparseTRMatrix<PetscVector>::printcontent(ConsoleOutput &output) const;
template<> void BlockGraphSparseTRMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const;

}
} /* end of namespace */

#endif
