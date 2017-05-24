#ifndef PASC_PETSCVECTOR_IMAGEDATA_H
#define	PASC_PETSCVECTOR_IMAGEDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/imagedata.h"
#include "external/petscvector/data/tsdata.h"

namespace pascinference {
namespace data {

template<> ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, std::string filename_data, int width, int height);
template<> void ImageData<PetscVector>::saveImage(std::string filename, bool save_original) const;



}
} /* end namespace */

#endif
