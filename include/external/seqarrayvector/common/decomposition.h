#ifndef PASC_SEQARRAY_COMMON_DECOMPOSITION_H
#define	PASC_SEQARRAY_COMMON_DECOMPOSITION_H

#include "external/seqarrayvector/algebra/vector/generalvector.h"
#include "general/common/decomposition.h"

namespace pascinference {
namespace algebra {

template<> Decomposition<SeqArrayVector>::Decomposition(int T, int R, int K, int xdim, int DDT_size);
template<> Decomposition<SeqArrayVector>::Decomposition(int T, BGMGraph<SeqArrayVector> &new_graph, int K, int xdim, int DDR_size);
template<> Decomposition<SeqArrayVector>::Decomposition(int T, BGMGraph<SeqArrayVector> &new_graph, int K, int xdim, int DDT_size, int DDR_size);

template<> Decomposition<SeqArrayVector>::~Decomposition();
template<> void Decomposition<SeqArrayVector>::compute_rank();
template<> void Decomposition<SeqArrayVector>::set_graph(BGMGraph<SeqArrayVector> &new_graph, int DDR_size);
/*
template<> void Decomposition<SeqArrayVector>::createGlobalVec_gamma(Vec *x_Vec) const;
template<> void Decomposition<SeqArrayVector>::createGlobalVec_data(Vec *x_Vec) const;
template<> void Decomposition<SeqArrayVector>::permute_TRblocksize(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const; 
template<> void Decomposition<SeqArrayVector>::permute_TRxdim(Vec orig_Vec, Vec new_Vec, bool invert) const;
*/

}
} /* end of namespace */

#endif
