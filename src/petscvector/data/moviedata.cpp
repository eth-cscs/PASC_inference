#include "external/petscvector/data/moviedata.h"

namespace pascinference {
namespace data {

template<>
MovieData<PetscVector>::MovieData(Decomposition<PetscVector> &new_decomposition, int width, int height, std::string filename_data, int type){
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
	/* type=0 -> TRn; type=1 -> TnR; type=2 -> nTR; */
	if(type == 0){
        this->decomposition->permute_TRb_to_dTRb(datapreload_Vec, data_Vec, decomposition->get_xdim(),false);
    }
	if(type == 1){
        this->decomposition->permute_TbR_to_dTRb(datapreload_Vec, data_Vec, decomposition->get_xdim(),false);
    }
	if(type == 2){
        this->decomposition->permute_bTR_to_dTRb(datapreload_Vec, data_Vec, decomposition->get_xdim(),false);
    }

	/* destroy preloaded vector */
//	TRYCXX( VecDestroy(&datapreload_Vec) );

	/* other vectors will be prepared after setting the model */
	this->destroy_datavector = true;
	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

template<>
MovieData<PetscVector>::MovieData(Decomposition<PetscVector> &new_decomposition, int width, int height){
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

template<>
void MovieData<PetscVector>::saveMovie_datavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	Timer timer_saveMovie;
	timer_saveMovie.restart();
	timer_saveMovie.start();

	std::ostringstream oss_name_of_file;

	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
	this->decomposition->createGlobalVec_data(&datasave_Vec);
	GeneralVector<PetscVector> datasave(datasave_Vec);

	/* save datavector - just for fun; to see if it was loaded in a right way */
	oss_name_of_file << "results/" << filename << "_datavector.bin";
	if(type == 0){
        this->decomposition->permute_TRb_to_dTRb(datasave_Vec, datavector->get_vector(), decomposition->get_xdim(), true);
    }
	if(type == 1){
        this->decomposition->permute_TbR_to_dTRb(datasave_Vec, datavector->get_vector(), decomposition->get_xdim(), true);
    }
	if(type == 2){
        this->decomposition->permute_bTR_to_dTRb(datasave_Vec, datavector->get_vector(), decomposition->get_xdim(), true);
    }

	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_saveMovie.stop();
	coutAll <<  " - Movie datavector saved in: " << timer_saveMovie.get_value_sum() << std::endl;
	coutAll.synchronize();

	LOG_FUNC_END
}

template<>
void MovieData<PetscVector>::saveMovie_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN

	Timer timer_saveMovie;
	timer_saveMovie.restart();
	timer_saveMovie.start();

	std::ostringstream oss_name_of_file;

	oss_name_of_file << "results/" << filename << "_gamma.bin";

	Vec gammasave_Vec;
    TRYCXX( VecDuplicate(gammavector->get_vector(), &gammasave_Vec) );
	this->decomposition->permute_TbR_to_dTRb(gammasave_Vec, gammavector->get_vector(), decomposition->get_K(), true);
	GeneralVector<PetscVector> gammasave(gammasave_Vec);
	gammasave.save_binary(oss_name_of_file.str());

	oss_name_of_file.str("");

	timer_saveMovie.stop();
	coutAll <<  " - Movie gammavector saved in: " << timer_saveMovie.get_value_sum() << std::endl;
	coutAll.synchronize();

	LOG_FUNC_END
}

template<>
void MovieData<PetscVector>::saveMovie_reconstructed(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	Timer timer_saveMovie;
	timer_saveMovie.restart();
	timer_saveMovie.start();

	std::ostringstream oss_name_of_file;

	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
	this->decomposition->createGlobalVec_data(&datasave_Vec);
	GeneralVector<PetscVector> datasave(datasave_Vec);

	/* compute recovered image */
	Vec gammak_Vec;
	IS gammak_is;
	Vec datan_recovered_Vec;
	IS datan_is;

	Vec data_recovered_Vec;
	TRYCXX( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRYCXX( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();
	int xdim = this->get_xdim();

	for(int k=0;k<K;k++){
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		for(int n=0;n<xdim;n++){
			/* get datan */
			this->decomposition->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_recovered_Vec, datan_is, &datan_recovered_Vec) );

			/* add to recovered image */
			TRYCXX( VecAXPY(datan_recovered_Vec, theta_arr[k*xdim + n], gammak_Vec) );

			/* restore data */
			TRYCXX( VecRestoreSubVector(data_recovered_Vec, datan_is, &datan_recovered_Vec) );
			TRYCXX( ISDestroy(&datan_is) );

		}

		TRYCXX( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";

	/* but at first, permute recovered data, datasave can be used */
	if(type == 0){
        this->decomposition->permute_TRb_to_dTRb(datasave_Vec, data_recovered_Vec, decomposition->get_xdim(), true);
    }
	if(type == 1){
        this->decomposition->permute_TbR_to_dTRb(datasave_Vec, data_recovered_Vec, decomposition->get_xdim(), true);
    }
	if(type == 2){
        this->decomposition->permute_bTR_to_dTRb(datasave_Vec, data_recovered_Vec, decomposition->get_xdim(), true);
    }

	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_saveMovie.stop();
	coutAll <<  " - Movie reconstructed saved in: " << timer_saveMovie.get_value_sum() << std::endl;
	coutAll.synchronize();

	LOG_FUNC_END
}


template<>
double MovieData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const {
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

	/* vector of absolute error */
	Vec data_abserr_Vec;
	TRYCXX( VecDuplicate(data_Vec, &data_abserr_Vec) );

	/* compute recovered signal */
	Vec datan_abserr_Vec;
	IS datan_is;
	Vec gammak_Vec;
	IS gammak_is;

	/* abserr = -solution */
	TRYCXX( VecCopy(solution_Vec,data_abserr_Vec));
	TRYCXX( VecScale(data_abserr_Vec,-1.0));

	double *theta_arr;
	TRYCXX( VecGetArray(this->thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();
	int xdim = this->get_xdim();

	for(int k=0;k<K;k++){
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );

		for(int n=0;n<xdim;n++){
			/* get datan */
			this->decomposition->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_abserr_Vec, datan_is, &datan_abserr_Vec) );

			/* add to recovered image */
			TRYCXX( VecAXPY(datan_abserr_Vec, theta_arr[k*xdim + n], gammak_Vec) );

			/* restore data */
			TRYCXX( VecRestoreSubVector(data_abserr_Vec, datan_is, &datan_abserr_Vec) );
			TRYCXX( ISDestroy(&datan_is) );

		}

		TRYCXX( VecRestoreSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	TRYCXX( VecRestoreArray(this->thetavector->get_vector(),&theta_arr) );

	/* compute mean(abs(solution - data_recovered) */
//	TRYCXX( VecAbs(data_abserr_Vec) );
//	TRYCXX( VecSum(data_abserr_Vec, &abserr) );
//	int T = this->get_T();
//	abserr = abserr/(double)T;

	TRYCXX( VecNorm(data_abserr_Vec, NORM_1, &abserr) );

	LOG_FUNC_END

	return abserr;
}



}
} /* end namespace */

