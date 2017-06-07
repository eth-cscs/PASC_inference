#ifndef PASC_PETSCVECTOR_SIGNALDATA_H
#define	PASC_PETSCVECTOR_SIGNALDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/signaldata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> SignalData<PetscVector>::SignalData(std::string filename_data);
template<> void SignalData<PetscVector>::set_decomposition(Decomposition<PetscVector> &new_decomposition);
template<> void SignalData<PetscVector>::saveSignal(std::string filename, bool save_original) const;
template<> double SignalData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const;

}
} /* end namespace */

#endif
