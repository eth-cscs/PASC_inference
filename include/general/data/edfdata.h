/** @file edfdata.h
 *  @brief 
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_EDFDATA_H
#define	PASC_EDFDATA_H

#include <iostream>
#include "general/common/common.h"
#include "general/model/tsmodel.h"
#include "general/data/tsdata.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class EdfData: public TSData<VectorBase> {
	protected:
		/* informations from EDF file */
		struct Record {
			std::string *hdr_label;
			std::string *hdr_transducer;
			std::string *hdr_units;
			double hdr_physicalMin; // TODO: double or int?
			double hdr_physicalMax;
			double hdr_digitalMin;
			double hdr_digitalMax;
			std::string *hdr_prefilter;
			int hdr_samples;
		};
		int hdr_ver;
		std::string hdr_patientID;
		std::string hdr_recordID;
		std::string hdr_startdate;
		std::string hdr_starttime;
		int hdr_bytes;
		int hdr_records;
		int hdr_duration;
		int hdr_ns;
		Record *hdr_records_detail;
		bool free_hdr_records_detail;

		void edfRead(std::string filename, int max_record_nmb = -1);

		/* preliminary data */
		int Tpreliminary;
		GeneralVector<VectorBase> *datavectorpreliminary;

	public:
		EdfData(std::string filename_data, int max_record_nmb = -1);
		~EdfData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveVTK(std::string filename) const;
		void saveVector(std::string filename, bool save_original) const;

		int get_Tpreliminary() const;
		void set_decomposition(Decomposition<VectorBase> &decomposition);

};


}
} /* end of namespace */

/* ------------- implementation ----------- */

namespace pascinference {
namespace data {

template<class VectorBase>
void EdfData<VectorBase>::edfRead(std::string filename, int max_record_nmb){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* set decomposition - from preliminary to real data */
template<class VectorBase>
void EdfData<VectorBase>::set_decomposition(Decomposition<VectorBase> &new_decomposition) {
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* from filename */
template<class VectorBase>
EdfData<VectorBase>::EdfData(std::string filename_data, int max_record_nmb){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EdfData<VectorBase>::~EdfData(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void EdfData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - version of this data format:            " << hdr_ver << std::endl;
	output <<  " - local patient identification:           " << hdr_patientID << std::endl;
	output <<  " - local recording identification:         " << hdr_recordID << std::endl;
	output <<  " - startdate of recording (dd.mm.yy):      " << hdr_startdate << std::endl;
	output <<  " - starttime of recording (hh.mm.ss):      " << hdr_starttime << std::endl;
	output <<  " - number of bytes in header record:       " << hdr_bytes << std::endl;
	output <<  " - number of data records (-1 if unknown): " << hdr_records << std::endl;
	output <<  " - duration of a data record, in seconds:  " << hdr_duration << std::endl;
	output <<  " - number of signals (ns) in data record:  " << hdr_ns << std::endl;
/*
	output <<  " - record details:" << std::endl;
	for(int i=0;i<hdr_ns;i++){
		output <<  "   - id:                 " << i << std::endl;
		output <<  "     label:              " << *hdr_records_detail[i].hdr_label << std::endl;
		output <<  "     transducer type:    " << *hdr_records_detail[i].hdr_transducer << std::endl;
		output <<  "     physical dimension: " << *hdr_records_detail[i].hdr_units << std::endl;
		output <<  "     physical minimum:   " << hdr_records_detail[i].hdr_physicalMin << std::endl;
		output <<  "     physical maximum:   " << hdr_records_detail[i].hdr_physicalMax << std::endl;
		output <<  "     digital minimum:    " << hdr_records_detail[i].hdr_digitalMin << std::endl;
		output <<  "     digital maximum:    " << hdr_records_detail[i].hdr_digitalMax << std::endl;
		output <<  "     prefiltering:       " << *hdr_records_detail[i].hdr_prefilter << std::endl;
		output <<  "     nr of samples:      " << hdr_records_detail[i].hdr_samples << std::endl;
	}
*/
	output <<  "----------------------------------------------------------------" << std::endl;

	output <<  " - Tpreliminary: " << this->get_T() << std::endl;

	if(this->decomposition){
		output <<  " - T           : " << this->get_T() << std::endl;
		output <<  " - xdim        : " << this->get_xdim() << std::endl;
		output <<  " - K           : " << this->get_K() << std::endl;
		output <<  " - R           : " << this->get_R() << std::endl;
	}

	if(this->tsmodel){
		output <<  " - model       : " << this->tsmodel->get_name() << std::endl;
	} else {
		output <<  " - model       : NO" << std::endl;
	}
	
	output <<  " - datavector  : ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<  " - gammavector : ";
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
void EdfData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - version of this data format:            " << hdr_ver << std::endl;
	output_global <<  " - local patient identification:           " << hdr_patientID << std::endl;
	output_global <<  " - local recording identification:         " << hdr_recordID << std::endl;
	output_global <<  " - startdate of recording (dd.mm.yy):      " << hdr_startdate << std::endl;
	output_global <<  " - starttime of recording (hh.mm.ss):      " << hdr_starttime << std::endl;
	output_global <<  " - number of bytes in header record:       " << hdr_bytes << std::endl;
	output_global <<  " - number of data records (-1 if unknown): " << hdr_records << std::endl;
	output_global <<  " - duration of a data record, in seconds:  " << hdr_duration << std::endl;
	output_global <<  " - number of signals (ns) in data record:  " << hdr_ns << std::endl;
/*
	output_global <<  " - record details:" << std::endl;
	for(int i=0;i<hdr_ns;i++){
		output_global <<  "   - id:                 " << i << std::endl;
		output_global <<  "     label:              " << *hdr_records_detail[i].hdr_label << std::endl;
		output_global <<  "     transducer type:    " << *hdr_records_detail[i].hdr_transducer << std::endl;
		output_global <<  "     physical dimension: " << *hdr_records_detail[i].hdr_units << std::endl;
		output_global <<  "     physical minimum:   " << hdr_records_detail[i].hdr_physicalMin << std::endl;
		output_global <<  "     physical maximum:   " << hdr_records_detail[i].hdr_physicalMax << std::endl;
		output_global <<  "     digital minimum:    " << hdr_records_detail[i].hdr_digitalMin << std::endl;
		output_global <<  "     digital maximum:    " << hdr_records_detail[i].hdr_digitalMax << std::endl;
		output_global <<  "     prefiltering:       " << *hdr_records_detail[i].hdr_prefilter << std::endl;
		output_global <<  "     nr of samples:      " << hdr_records_detail[i].hdr_samples << std::endl;
	}
*/
	output_global <<  "----------------------------------------------------------------" << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - Tpreliminary: " << this->get_T() << std::endl;

	if(this->decomposition){
		output_global <<  " - T           : " << this->get_T() << std::endl;
		output_local  <<  "  - Tlocal     : " << this->get_Tlocal() << std::endl;
		output_local.synchronize();
		output_global <<  " - xdim        : " << this->get_xdim() << std::endl;
		output_global <<  " - K           : " << this->get_K() << std::endl;
		output_global <<  " - R           : " << this->get_R() << std::endl;
		output_local  <<  "  - Rlocal     : " << this->get_Rlocal() << std::endl;
		output_local.synchronize();
	}

	if(this->tsmodel){
		output_global <<  " - model       : " << this->tsmodel->get_name() << std::endl;
	} else {
		output_global <<  " - model       : NO" << std::endl;
	}
	
	output_global <<  " - datavector  : ";
	if(this->datavector){
		output_global << "YES (size: " << this->datavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->datavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global <<  " - gammavector : ";
	if(this->gammavector){
		output_global << "YES (size: " << this->gammavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->gammavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global << " - thetavector : ";
	if(this->thetavector){
		output_global << "YES (size: " << this->thetavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->thetavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void EdfData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print the content of the data */
	output <<  " - datavector  : ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gammavector : ";
	if(this->gammavector){
		output << *this->gammavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - thetavector : ";
	if(this->thetavector){
		output << *this->thetavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void EdfData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output_local <<  " - datavector : ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - gammavector : ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - thetavector : ";
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
std::string EdfData<VectorBase>::get_name() const {
	return "EDF Time-series Data";
}

template<class VectorBase>
int EdfData<VectorBase>::get_Tpreliminary() const{
	return this->Tpreliminary;
}

template<class VectorBase>
void EdfData<VectorBase>::saveVTK(std::string filename) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void EdfData<VectorBase>::saveVector(std::string filename, bool save_original) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


}
} /* end namespace */

#endif
