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

	//TODO
	
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

	//TODO
	
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

template<class VectorBase>
void Signal1DData<VectorBase>::saveSignal1D(std::string filename, bool save_original) const{
	LOG_FUNC_BEGIN
	
	//TODO
	
	LOG_FUNC_END
}

template<class VectorBase>
int Signal1DData<VectorBase>::get_Tpreliminary() const{
	return this->Tpreliminary;
}

template<class VectorBase>
double Signal1DData<VectorBase>::compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const {
	LOG_FUNC_BEGIN	

	//TODO
	
	LOG_FUNC_END
	
	return -1.0;
}



}
} /* end namespace */

#endif
