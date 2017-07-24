#include "external/petscvector/data/imagedata.h"

namespace pascinference {
namespace data {

template<>
ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, int width, int height, std::string filename_data){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;
	this->decomposition = &new_decomposition;

	/* prepare preliminary datavector and load data */
	Vec datapreload_Vec;
	this->decomposition->createGlobalVec_data(&datapreload_Vec);
	GeneralVector<PetscVector> datapreload(datapreload_Vec);
	datapreload.load_global(filename_data);

	/* prepare real datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);

	/* permute orig to new using parallel layout */
	this->decomposition->permute_gTbR_to_pdTRb(datapreload_Vec, data_Vec, decomposition->get_xdim(),false);

	/* destroy preloaded vector */
//	TRYCXX( VecDestroy(&datapreload_Vec) );

	/* other vectors will be prepared after setting the model */
	this->destroy_datavector = true;
	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

template<>
ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, int width, int height){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;
	this->decomposition = &new_decomposition;

	/* prepare datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);

	/* other vectors will be prepared after setting the model */
	this->destroy_datavector = true;
	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

}
} /* end namespace */

