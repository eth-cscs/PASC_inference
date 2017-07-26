#ifndef PASC_PETSCVECTOR_BLOCKGRAPHSPARSEMATRIX_H
#define	PASC_PETSCVECTOR_BLOCKGRAPHSPARSEMATRIX_H

#include "general/algebra/matrix/blockgraphsparse.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/matrix/generalmatrix.h"
#include "external/petscvector/common/decomposition.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BlockGraphSparseMatrix<PetscVector>::ExternalContent : public GeneralMatrix<PetscVector>::ExternalContent {
    public:
		virtual std::string get_name() const {
            return "BlockGraphSparseMatrix PETSc content";
		}
};

template<> BlockGraphSparseMatrix<PetscVector>::BlockGraphSparseMatrix(Decomposition<PetscVector> &new_decomposition, double alpha, GeneralVector<PetscVector> *new_coeffs);
template<> BlockGraphSparseMatrix<PetscVector>::~BlockGraphSparseMatrix();
template<> void BlockGraphSparseMatrix<PetscVector>::printcontent(ConsoleOutput &output) const;
template<> void BlockGraphSparseMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const;

template<> BlockGraphSparseMatrix<PetscVector>::ExternalContent * BlockGraphSparseMatrix<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
