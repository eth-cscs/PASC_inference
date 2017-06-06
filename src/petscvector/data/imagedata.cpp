#include "external/petscvector/data/imagedata.h"

namespace pascinference {
namespace data {

template<>
ImageData<PetscVector>::ImageData(Decomposition<PetscVector> &new_decomposition, std::string filename_data, int width, int height){
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
	this->decomposition->permute_TRxdim(datapreload_Vec, data_Vec);
	
	/* destroy preloaded vector */
//	TRYCXX( VecDestroy(&datapreload_Vec) );
	
	/* other vectors will be prepared after setting the model */
	this->destroy_datavector = true;
	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

template<>
void ImageData<PetscVector>::saveImage(std::string filename, bool save_original) const{
	LOG_FUNC_BEGIN
	
	Timer timer_saveImage; 
	timer_saveImage.restart();
	timer_saveImage.start();

	std::ostringstream oss_name_of_file;

	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
	this->decomposition->createGlobalVec_data(&datasave_Vec);
	GeneralVector<PetscVector> datasave(datasave_Vec);

	Vec gammasave_Vec;
	this->decomposition->createGlobalVec_gamma(&gammasave_Vec);
	GeneralVector<PetscVector> gammasave(gammasave_Vec);

	/* save datavector - just for fun; to see if it was loaded in a right way */
	if(save_original){
		oss_name_of_file << "results/" << filename << "_original.bin";
		this->decomposition->permute_TRxdim(datasave_Vec, datavector->get_vector(), true);
		datasave.save_binary(oss_name_of_file.str());
		oss_name_of_file.str("");
	}

	/* save gamma */
	oss_name_of_file << "results/" << filename << "_gamma.bin";
	this->decomposition->permute_TRK(gammasave_Vec, gammavector->get_vector(), true);
	gammasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* compute recovered image */
	Vec gammak_Vec;
	IS gammak_is;

	Vec data_recovered_Vec;
	TRYCXX( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRYCXX( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();

	for(int k=0;k<K;k++){ 
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		TRYCXX( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		TRYCXX( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}	

	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";
	
	/* but at first, permute recovered data, datasave can be used */
	this->decomposition->permute_TRxdim(datasave_Vec, data_recovered_Vec, true);
	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* destroy vectors with original layout */
//	TRYCXX( VecDestroy(&datasave_Vec) );

	timer_saveImage.stop();
	coutAll <<  " - problem saved in: " << timer_saveImage.get_value_sum() << std::endl;
	coutAll.synchronize();
	
	LOG_FUNC_END
}

template<>
double ImageData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const {
	LOG_FUNC_BEGIN	

	double abserr;

	Vec data_Vec = this->datavector->get_vector();
	Vec solution_Vec = solution.get_vector();
	Vec gammavector_Vec = this->gammavector->get_vector();

	/* transfer data to GPU */
	#ifdef USE_CUDA
		cuda_copytoGPU(data_Vec);
		cuda_copytoGPU(solution_Vec);
		cuda_copytoGPU(gammavector_Vec);
	#endif

	/* compute recovered signal */
	Vec gammak_Vec;
	IS gammak_is;

	Vec data_abserr_Vec;
	TRYCXX( VecDuplicate(data_Vec, &data_abserr_Vec) );
	
	/* abserr = -solution */
	TRYCXX( VecCopy(solution_Vec,data_abserr_Vec)); 
	TRYCXX( VecScale(data_abserr_Vec,-1.0));

	double *theta_arr;
	TRYCXX( VecGetArray(this->thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();

	for(int k=0;k<K;k++){ 
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );

		/* add to recovered image */
		TRYCXX( VecAXPY(data_abserr_Vec, theta_arr[k], gammak_Vec) );

		TRYCXX( VecRestoreSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}	

	TRYCXX( VecRestoreArray(this->thetavector->get_vector(),&theta_arr) );

	/* compute mean(abs(solution - data_recovered) */
//	TRYCXX( VecAbs(data_abserr_Vec) );
//	TRYCXX( VecSum(data_abserr_Vec, &abserr) );
//	int T = this->get_T();
//	abserr = abserr/(double)T;

	TRYCXX( VecNorm(data_abserr_Vec, NORM_2, &abserr) );

	LOG_FUNC_END
	
	return abserr;
}



}
} /* end namespace */

