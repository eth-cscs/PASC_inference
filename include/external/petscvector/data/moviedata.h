#ifndef PASC_PETSCVECTOR_MOVIEDATA_H
#define	PASC_PETSCVECTOR_MOVIEDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/moviedata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> MovieData<PetscVector>::MovieData(Decomposition<PetscVector> &new_decomposition, std::string filename_data, int width, int height);
template<> void MovieData<PetscVector>::saveMovie(std::string filename, bool save_original) const;
template<> double MovieData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const;

}
} /* end namespace */

#endif
