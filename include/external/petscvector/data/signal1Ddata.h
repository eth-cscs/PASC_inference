#ifndef PASC_PETSCVECTOR_SIGNAL1DDATA_H
#define	PASC_PETSCVECTOR_SIGNAL1DDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/signal1Ddata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> Signal1DData<PetscVector>::Signal1DData(std::string filename_data);
template<> void Signal1DData<PetscVector>::set_decomposition(Decomposition<PetscVector> &new_decomposition);
template<> void Signal1DData<PetscVector>::saveSignal1D(std::string filename, bool save_original) const;
template<> double Signal1DData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const;

}
} /* end namespace */

#endif
