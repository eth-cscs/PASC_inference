#ifndef PASC_PETSCVECTOR_IMAGEDATA_H
#define	PASC_PETSCVECTOR_IMAGEDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/imagedata.h"
#include "external/petscvector/data/tsdata.h"
#include "external/petscvector/common/common.h"

namespace pascinference {
namespace data {

template<> ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, std::string filename_data, int width, int height);
template<> ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, int width, int height);

template<> void ImageData<PetscVector>::saveImage_datavector(std::string filename) const;
template<> void ImageData<PetscVector>::saveImage_gammavector(std::string filename) const;
template<> void ImageData<PetscVector>::saveImage_reconstructed(std::string filename) const;

template<> double ImageData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const;

}
} /* end namespace */

#endif
