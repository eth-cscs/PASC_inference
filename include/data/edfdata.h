/** @file edfdata.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_EDFDATA_H
#define	PASC_EDFDATA_H

#ifndef USE_PETSCVECTOR
 #error 'EDFDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "data/tsdata.h"
#include "model/tsmodel.h"

namespace pascinference {

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

		int R;
		int T;
		int Tlocal;
		
		void edfRead(std::string filename);
	public:
		EdfData(std::string filename_data);
		~EdfData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		int get_R() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<>
void EdfData<PetscVector>::edfRead(std::string filename){
	LOG_FUNC_BEGIN

	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	myfile.seekg(0, std::ios::beg);

	int i;
	char buffer[100];

	/* ------ HEADER ------ */
	myfile.read(buffer, 8);
	hdr_ver = atoi(buffer);

	myfile.read(buffer, 80);
	hdr_patientID = std::string(buffer);

	myfile.read(buffer, 80);
	hdr_recordID = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_startdate = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_starttime = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_bytes = atoi(buffer);

	myfile.read(buffer, 44);
	/* reserved */
	
	myfile.read(buffer, 8);
	hdr_records = atoi(buffer);

	myfile.read(buffer, 8);
	hdr_duration = atoi(buffer);

	myfile.read(buffer, 4);
	hdr_ns = atoi(buffer);

	/* arrays */
	hdr_records_detail = (Record*)malloc(sizeof(Record)*hdr_ns);
	free_hdr_records_detail = true;
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 16);
		hdr_records_detail[i].hdr_label = new std::string(buffer);
	}
			
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 80);
		hdr_records_detail[i].hdr_transducer = new std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_units = new std::string(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_physicalMin = atof(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_physicalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_digitalMin = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_digitalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 80);
		hdr_records_detail[i].hdr_prefilter = new std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_samples = atoi(buffer);
	}	

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 32);
		/* reserved */
	}	

	/* ------ PREPARE DATAVECTOR ------ */
	/* compute vector lengths */
	T = hdr_records_detail[0].hdr_samples*hdr_records;
	R = hdr_ns-1;

	/* for data prepare vector of length T and spit it into processors */
	Vec layout;

	/* try to make a global vector of length T and then get size of local portion */
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,this->T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&(this->Tlocal)) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */

	Vec datavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
	TRY( VecSetSizes(datavector_Vec,this->Tlocal*R,T*R) );
	TRY( VecSetFromOptions(datavector_Vec) );
	this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
	this->destroy_datavector = true;

	/* get ownership range */
	int t_begin, t_end, t_length; 
	TRY( VecGetOwnershipRange(datavector_Vec, &t_begin, &t_end) );
	t_begin = ((double)t_begin)/((double)R);
	t_end = ((double)t_end)/((double)R);			
	t_length = t_end - t_begin;

	/* ------ RECORDS ------ */
	int recnum, ii, samplei, index;
	double scalefac, dc;
    
	int16_t value;

    for(recnum = 0; recnum < hdr_records; recnum++){
		for(ii = 0; ii < R; ii++){
			scalefac = (hdr_records_detail[ii].hdr_physicalMax - hdr_records_detail[ii].hdr_physicalMin)/(double)(hdr_records_detail[ii].hdr_digitalMax - hdr_records_detail[ii].hdr_digitalMin);
			dc = hdr_records_detail[ii].hdr_physicalMax - scalefac*hdr_records_detail[ii].hdr_digitalMax;

			for(samplei=0; samplei < hdr_records_detail[ii].hdr_samples; samplei++){
				myfile.read((char *)&value, sizeof(int16_t)); /* read block of memory */
				value =  value * scalefac + dc;
				index = ii*this->T + recnum*hdr_records_detail[ii].hdr_samples + samplei;
				TRY( VecSetValue(datavector_Vec, index, value, INSERT_VALUES) );

			}
        }
    }

	/* vector is prepared */
	TRY( VecAssemblyBegin(datavector_Vec) );
	TRY( VecAssemblyEnd(datavector_Vec) );

	/* close file */
    myfile.close();		

	LOG_FUNC_END
}

/* from filename */
template<class VectorBase>
EdfData<VectorBase>::EdfData(std::string filename_data){
	LOG_FUNC_BEGIN

	/* read data from input file */
	edfRead(filename_data);

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EdfData<VectorBase>::~EdfData(){
	LOG_FUNC_BEGIN
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->free_hdr_records_detail){
		free(this->hdr_records_detail);
	}

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
void EdfData<VectorBase>::printcontent(ConsoleOutput &output) const {
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
void EdfData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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
std::string EdfData<VectorBase>::get_name() const {
	return "EDF Time-series Data";
}

template<class VectorBase>
int EdfData<VectorBase>::get_R() const{
	return this->R;
}



} /* end namespace */

#endif
