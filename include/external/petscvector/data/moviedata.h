#ifndef PASC_PETSCVECTOR_MOVIEDATA_H
#define	PASC_PETSCVECTOR_MOVIEDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/moviedata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> MovieData<PetscVector>::MovieData(Decomposition<PetscVector> &new_decomposition, int width, int height, std::string filename_data, int type);
template<> MovieData<PetscVector>::MovieData(Decomposition<PetscVector> &new_decomposition, int width, int height);

}
} /* end namespace */

#endif
