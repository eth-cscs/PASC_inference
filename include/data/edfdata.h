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
#include "data/tsdata_global.h"
#include "model/tsmodel_global.h"

namespace pascinference {

class EdfData: public TSData_Global {
	protected:
		/* informations from EDF file */
		struct Record {
			std::string hdr_label;
			std::string hdr_transducer;
			std::string hdr_units;
			double hdr_physicalMin; // TODO: double or int?
			double hdr_physicalMax;
			double hdr_digitalMin;
			double hdr_digitalMax;
			std::string hdr_prefilter;
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
		
		void read_from_file(std::ifstream &myfile, char *outchar, int length){
			myfile.read(outchar, length);
		};

		void edfRead(std::string filename);
	public:
//		EdfData();
//		EdfData(TSModel_Global &tsmodel) {};
		EdfData(std::string filename);
		~EdfData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/*EdfData::EdfData(){
	LOG_FUNC_BEGIN

	this->tsmodel = NULL;

	this->datavector = NULL;
	destroy_datavector = false;

	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	free_hdr_records_detail = false;
	
	LOG_FUNC_END
}
*/

void EdfData::edfRead(std::string filename){
	LOG_FUNC_BEGIN

	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	myfile.seekg(0, std::ios::beg);

	int i;
	char buffer[100];

	/* ------ HEADER ------ */
	read_from_file(myfile, buffer, 8);
	hdr_ver = atoi(buffer);

	read_from_file(myfile, buffer, 80);
	hdr_patientID = std::string(buffer);

	read_from_file(myfile, buffer, 80);
	hdr_recordID = std::string(buffer);

	read_from_file(myfile, buffer, 8);
	hdr_startdate = std::string(buffer);

	read_from_file(myfile, buffer, 8);
	hdr_starttime = std::string(buffer);

	read_from_file(myfile, buffer, 8);
	hdr_bytes = atoi(buffer);

	read_from_file(myfile, buffer, 44);
	/* reserved */
	
	read_from_file(myfile, buffer, 8);
	hdr_records = atoi(buffer);

	read_from_file(myfile, buffer, 8);
	hdr_duration = atoi(buffer);

	read_from_file(myfile, buffer, 4);
	hdr_ns = atoi(buffer);

	/* arrays */
	hdr_records_detail = (Record*)malloc(sizeof(Record)*hdr_ns);
	free_hdr_records_detail = true;

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 16);
		hdr_records_detail[i].hdr_label = std::string(buffer);
	}
			
	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 80);
		hdr_records_detail[i].hdr_transducer = std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_units = std::string(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_physicalMin = atof(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_physicalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_digitalMin = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_digitalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 80);
		hdr_records_detail[i].hdr_prefilter = std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 8);
		hdr_records_detail[i].hdr_samples = atoi(buffer);
	}	

	for(i=0;i<hdr_ns;i++){
		read_from_file(myfile, buffer, 32);
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

//    read_from_file(myfile, &value, sizeof(uint16_t));
//	value = atoi(buffer);

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

/*    
    record = zeros(hdr.ns, hdr.samples(1)*hdr.records);
    
    for ii = 1:numel(hdr.label)
        ctr = 1;
        for jj = 1:hdr.records
            try
                record(ii, ctr : ctr + hdr.samples(ii) - 1) = tmpdata(jj).data{ii};
            end
            ctr = ctr + length(tmpdata(jj).data{ii});
        end
    end
*/

//	coutMaster << "size of int16:" << sizeof(int16_t) << std::endl;

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	/* close file */
    myfile.close();		

	LOG_FUNC_END
}

/* from filename */
EdfData::EdfData(std::string filename){
	LOG_FUNC_BEGIN

	/* read data from input file */
	edfRead(filename);

	LOG_FUNC_END
}

/* destructor */
EdfData::~EdfData(){
	LOG_FUNC_BEGIN
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->free_hdr_records_detail){
		free(this->hdr_records_detail);
	}

	LOG_FUNC_END
}


/* print info about data */
void EdfData::print(ConsoleOutput &output) const {
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

	output <<  " - record details:" << std::endl;
	for(int i=0;i<hdr_ns;i++){
		output <<  "   - id:                 " << i << std::endl;
		output <<  "     label:              " << hdr_records_detail[i].hdr_label << std::endl;
		output <<  "     transducer type:    " << hdr_records_detail[i].hdr_transducer << std::endl;
		output <<  "     physical dimension: " << hdr_records_detail[i].hdr_units << std::endl;
		output <<  "     physical minimum:   " << hdr_records_detail[i].hdr_physicalMin << std::endl;
		output <<  "     physical maximum:   " << hdr_records_detail[i].hdr_physicalMax << std::endl;
		output <<  "     digital minimum:    " << hdr_records_detail[i].hdr_digitalMin << std::endl;
		output <<  "     digital maximum:    " << hdr_records_detail[i].hdr_digitalMax << std::endl;
		output <<  "     prefiltering:       " << hdr_records_detail[i].hdr_prefilter << std::endl;
		output <<  "     nr of samples:      " << hdr_records_detail[i].hdr_samples << std::endl;
	}

	output <<  "----------------------------------------------------------------" << std::endl;

/*	output <<  " - T:           " << this->get_T() << std::endl;
	output <<  " - xdim:        " << this->get_xdim() << std::endl;
	output <<  " - K:           " << this->tsmodel->get_K() << std::endl;
	output <<  " - model:       " << this->tsmodel->get_name() << std::endl;
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
*/
	output.synchronize();

	LOG_FUNC_END
}

/* print info about data */
void EdfData::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN
	
	/* give information about presence of the data */
	output_global <<  " - T:           " << this->get_T() << std::endl;
	output_local  <<  "  - Tlocal:     " << this->tsmodel->get_Tlocal() << std::endl;
	output_local.synchronize();

	output_global <<  " - xdim:        " << this->get_xdim() << std::endl;
	output_global <<  " - K:           " << this->tsmodel->get_K() << std::endl;

	output_global <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	
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
void EdfData::printcontent(ConsoleOutput &output) const {
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
void EdfData::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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

std::string EdfData::get_name() const {
	return "EDF Time-series Data";
}





} /* end namespace */

#endif
