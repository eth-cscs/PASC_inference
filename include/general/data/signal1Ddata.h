/** @file signal1Ddata.h
 *  @brief class for manipulation with one dimensional data
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIGNAL1DDATA_H
#define	PASC_SIGNAL1DDATA_H

#include <iostream>
#include "general/common/common.h"
#include "general/model/tsmodel.h"
#include "general/data/tsdata.h"

namespace pascinference {
namespace data {

/** class Signal1DData
 * @brief data of one-dimensional signal
 * 
 * Class for manipulation with data from simple one-dimensional signal.
 */ 
template<class VectorBase>
class Signal1DData: public TSData<VectorBase> {
	protected:
		/* preliminary data */
		int Tpreliminary;
		GeneralVector<VectorBase> *datavectorpreliminary;

	public:
		Signal1DData(std::string filename_data);
		~Signal1DData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveSignal1D(std::string filename, bool save_original=true) const;

		int get_Tpreliminary() const;
		void set_decomposition(Decomposition<VectorBase> &decomposition);
		double compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace data {

/* from filename */
template<class VectorBase>
Signal1DData<VectorBase>::Signal1DData(std::string filename_data){
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

/* destructor */
template<class VectorBase>
Signal1DData<VectorBase>::~Signal1DData(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END
}


/* set decomposition - from preliminary to real data */
template<class VectorBase>
void Signal1DData<VectorBase>::set_decomposition(Decomposition<VectorBase> &new_decomposition) {
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;

	/* prepare real datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	this->destroy_datavector = true;
	
	/* permute orig to new using parallel layout */
	Vec datapreload_Vec = datavectorpreliminary->get_vector();
	this->decomposition->permute_TRxdim(datapreload_Vec, data_Vec);
	
	/* destroy preliminary data */
	TRYCXX(VecDestroy(&datapreload_Vec));
	
	LOG_FUNC_END
}

/* print info about data */
template<class VectorBase>
void Signal1DData<VectorBase>::print(ConsoleOutput &output) const {
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
void Signal1DData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
void Signal1DData<VectorBase>::printcontent(ConsoleOutput &output) const {
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
void Signal1DData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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
std::string Signal1DData<VectorBase>::get_name() const {
	return "Signal1D Time-series Data";
}

template<>
void Signal1DData<PetscVector>::saveSignal1D(std::string filename, bool save_original) const{
	Timer timer_saveSignal1D; 
	timer_saveSignal1D.restart();
	timer_saveSignal1D.start();

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

	/* compute recovered signal */
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

	timer_saveSignal1D.stop();
	coutAll <<  " - problem saved in: " << timer_saveSignal1D.get_value_sum() << std::endl;
	coutAll.synchronize();
}

template<class VectorBase>
int Signal1DData<VectorBase>::get_Tpreliminary() const{
	return this->Tpreliminary;
}

template<class VectorBase>
double Signal1DData<VectorBase>::compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const {
	LOG_FUNC_BEGIN	

	double abserr;

	Vec data_Vec = this->datavector->get_vector();
	Vec solution_Vec = solution.get_vector();
	Vec gammavector_Vec = this->gammavector->get_vector();

	/* transfer data to GPU */
	#ifdef USE_CUDA
		TRYCXX( VecCUDACopyToGPU(data_Vec) );
		TRYCXX( VecCUDACopyToGPU(solution_Vec) );
		TRYCXX( VecCUDACopyToGPU(gammavector_Vec) );
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

	TRYCXX( VecNorm(data_abserr_Vec, NORM_1, &abserr) );
	
	LOG_FUNC_END
	
	return abserr;
}



}
} /* end namespace */

#endif