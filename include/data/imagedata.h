/** @file imagedata.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_IMAGEDATA_H
#define	PASC_IMAGEDATA_H

#ifndef USE_PETSCVECTOR
 #error 'IMAGEDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "matrix/blockgraph.h"
#include "model/tsmodel.h"
#include "data/tsdata.h"

namespace pascinference {

template<class VectorBase>
class ImageData: public TSData<VectorBase> {
	protected:
		int R;
	public:
		ImageData(std::string filename_data);
		~ImageData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveImage(std::string filename) const;

		int get_R() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* from filename */
template<class VectorBase>
ImageData<VectorBase>::ImageData(std::string filename_data){
	LOG_FUNC_BEGIN

	/* ------ PREPARE DATAVECTOR ------ */
	Vec datavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
	TRY( VecSetFromOptions(datavector_Vec) );

	this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
	this->destroy_datavector = true;

	/* load image from file */
	this->datavector->load_local(filename_data);

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	/* compute vector lengths */
	this->T = 1;
	TRY(VecGetSize(datavector_Vec,&(this->R)));

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
		output <<  " - K:           " << this->tsmodel->get_K() << std::endl;
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
		output_local  <<  "  - Tlocal:     " << this->tsmodel->get_Tlocal() << std::endl;
		output_local.synchronize();

		output_global <<  " - xdim:        " << this->get_xdim() << std::endl;
		output_global <<  " - K:           " << this->tsmodel->get_K() << std::endl;

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

template<class VectorBase>
int ImageData<VectorBase>::get_R() const{
	return this->R;
}

template<>
void ImageData<PetscVector>::saveImage(std::string filename) const{
	Timer timer_saveImage; 
	timer_saveImage.restart();
	timer_saveImage.start();

	std::ostringstream oss_name_of_file;

	/* save datavector - just for fun; to see if it was loaded in a right way */
	oss_name_of_file << "results/" << filename << "_original.bin";
	datavector->save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* save gammas and compute recovered image */
	Vec gamma_Vec = gammavector->get_vector();
	Vec gammak_Vec;
	GeneralVector<PetscVector> *gammak;
	IS gammak_is;

	Vec data_recovered_Vec;
	TRY( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRY( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRY( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = get_K();
	int Tlocal = get_Tlocal();
	int k;
	for(k=0;k<K;k++){ 
		/* get gammak */
		TRY( ISCreateStride(PETSC_COMM_WORLD, R*Tlocal, Tbegin*K*R + k*Tlocal*R, 1, &gammak_is) );
		TRY( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		gammak = new GeneralVector<PetscVector>(gammak_Vec);

		/* save gammak */
		oss_name_of_file << "results/" << filename << "_gamma" << k << ".bin";
		gammak->save_binary(oss_name_of_file.str());
		oss_name_of_file.str("");

		oss_name_of_file << "results/" << filename << "_gamma" << k << ".txt";
		gammak->save_ascii(oss_name_of_file.str());
		oss_name_of_file.str("");


		/* add to recovered image */
		TRY( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		free(gammak);
	
		TRY( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		TRY( ISDestroy(&gammak_is) );
	}	

	TRY( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";
	data_recovered.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_saveImage.stop();
	coutAll <<  " - problem saved in: " << timer_saveImage.get_value_sum() << std::endl;
	coutAll.synchronize();
}


} /* end namespace */

#endif
