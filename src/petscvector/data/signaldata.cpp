#include "external/petscvector/data/signaldata.h"

namespace pascinference {
namespace data {

/* from filename */
template<>
SignalData<PetscVector>::SignalData(std::string filename_data){
	LOG_FUNC_BEGIN

	/* prepare preliminary datavector and load data */
	Vec datapreload_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&datapreload_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(datapreload_Vec, VECMPICUDA));
	#endif

	this->datavectorpreliminary = new GeneralVector<PetscVector>(datapreload_Vec);
	this->datavectorpreliminary->load_global(filename_data);

	/* get the size of the loaded vector */
	TRYCXX( VecGetSize(datapreload_Vec, &Tpreliminary) );

	/* other vectors will be prepared after setting the model */
	this->destroy_datavector = true;
	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

template<>
void SignalData<PetscVector>::set_decomposition(Decomposition<PetscVector> &new_decomposition, int type) {
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;

	/* prepare real datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	this->destroy_datavector = true;

	/* permute orig to new using parallel layout */
	Vec datapreload_Vec = datavectorpreliminary->get_vector();
	this->decomposition->permute_to_pdTRb(datapreload_Vec, data_Vec, decomposition->get_xdim(), type, false);

	/* destroy preliminary data */
	TRYCXX(VecDestroy(&datapreload_Vec));

	LOG_FUNC_END
}

}
} /* end namespace */

