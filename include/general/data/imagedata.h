/** @file imagedata.h
 *  @brief class for manipulation with image data
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_IMAGEDATA_H
#define	PASC_IMAGEDATA_H

#include <iostream>
#include "general/common/common.h"
#include "general/algebra/graph/bgmgraph.h"
#include "general/model/tsmodel.h"
#include "general/data/tsdata.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class ImageData: public TSData<VectorBase> {
	protected:
		int width;
		int height;
	public:
		ImageData(Decomposition<VectorBase> &decomposition, std::string filename_data, int width, int height);
		~ImageData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveImage(std::string filename, bool save_original=true) const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace data {

/* from filename */
template<class VectorBase>
ImageData<VectorBase>::ImageData(Decomposition<VectorBase> &new_decomposition, std::string filename_data, int width, int height){
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

/* destructor */
template<class VectorBase>
ImageData<VectorBase>::~ImageData(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void ImageData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	if(this->tsmodel){
		output <<  " - T:           " << this->get_T() << std::endl;
		output <<  " - xdim:        " << this->get_xdim() << std::endl;
		output <<  " - K:           " << this->get_K() << std::endl;
		output <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	} else {
		output <<  " - model:       NO" << std::endl;
	}
	output <<  " - R:           " << this->get_R() << std::endl;
	
	output <<  " - datavector:  ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - gammavector: ";
	if(this->gammavector){
		output << "YES (size: " << this->gammavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - thetavector: ";
	if(this->thetavector){
		output << "YES (size: " << this->thetavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}

	output.synchronize();

	LOG_FUNC_END
}

/* print info about data */
template<class VectorBase>
void ImageData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	if(this->tsmodel){
		output_global <<  " - T:           " << this->get_T() << std::endl;
		output_local  <<  "  - Tlocal:     " << this->get_Tlocal() << std::endl;
		output_local.synchronize();

		output_global <<  " - xdim:        " << this->get_xdim() << std::endl;
		output_global <<  " - K:           " << this->get_K() << std::endl;

		output_global <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	} else {
		output_global <<  " - model:       NO" << std::endl;
	}
	output_global <<  " - R:           " << this->get_R() << std::endl;
	
	output_global <<  " - datavector:  ";
	if(this->datavector){
		output_global << "YES (size: " << this->datavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->datavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global <<   " - gammavector: ";
	if(this->gammavector){
		output_global << "YES (size: " << this->gammavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->gammavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global <<   " - thetavector: ";
	if(this->thetavector){
		output_global << "YES (size: " << this->thetavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->thetavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void ImageData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output <<  " - datavector: ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gammavector: ";
	if(this->gammavector){
		output << *this->gammavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - thetavector: ";
	if(this->thetavector){
		output << *this->thetavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void ImageData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output_local <<  " - datavector: ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - gammavector: ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - thetavector: ";
	if(this->thetavector){
		output_local << *this->thetavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
std::string ImageData<VectorBase>::get_name() const {
	return "Image Time-series Data";
}

template<>
void ImageData<PetscVector>::saveImage(std::string filename, bool save_original) const{
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
}


}
} /* end namespace */

#endif