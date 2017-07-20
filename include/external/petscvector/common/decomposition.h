#ifndef PASC_PETSCVECTOR_COMMON_DECOMPOSITION_H
#define	PASC_PETSCVECTOR_COMMON_DECOMPOSITION_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"
#include "general/common/decomposition.h"

namespace pascinference {
namespace algebra {

template<> Decomposition<PetscVector>::Decomposition(int T, int R, int K, int xdim, int DDT_size);
template<> Decomposition<PetscVector>::Decomposition(int T, BGMGraph<PetscVector> &new_graph, int K, int xdim, int DDR_size);
template<> Decomposition<PetscVector>::Decomposition(int T, BGMGraph<PetscVector> &new_graph, int K, int xdim, int DDT_size, int DDR_size);

template<> Decomposition<PetscVector>::~Decomposition();
template<> void Decomposition<PetscVector>::compute_rank();
template<> void Decomposition<PetscVector>::set_graph(BGMGraph<PetscVector> &new_graph, int DDR_size);
template<> void Decomposition<PetscVector>::set_new_graph(BGMGraph<PetscVector> &new_graph, int DDR_size);

template<> void Decomposition<PetscVector>::createGlobalVec_gamma(Vec *x_Vec) const;
template<> void Decomposition<PetscVector>::createGlobalVec_data(Vec *x_Vec) const;

template<> void Decomposition<PetscVector>::permute_to_pdTRb(Vec orig_Vec, Vec new_Vec, int blocksize, int type, bool invert) const;
template<> void Decomposition<PetscVector>::permute_gTRb_to_pdTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const;
template<> void Decomposition<PetscVector>::permute_gTbR_to_pdTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const;
template<> void Decomposition<PetscVector>::permute_gbTR_to_pdTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const;

/* PETSc specific stuff */

template<> void Decomposition<PetscVector>::createIS_dTR_to_pdTRb(IS *is, int blocksize) const;
template<> void Decomposition<PetscVector>::createIS_gammaK(IS *is, int k) const;
template<> void Decomposition<PetscVector>::createIS_datan(IS *is, int n) const;

}
} /* end of namespace */

#endif
