#ifndef PASC_PETSCVECTOR_IMAGEDATA_H
#define	PASC_PETSCVECTOR_IMAGEDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/imagedata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, int width, int height, std::string filename_data);
template<> ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, int width, int height);

}
} /* end namespace */

#endif
