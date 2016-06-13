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

	LOG_FUNC_END
}
*/

/* from filename */
EdfData::EdfData(std::string filename){
	LOG_FUNC_BEGIN


	LOG_FUNC_END
}

/* destructor */
EdfData::~EdfData(){
	LOG_FUNC_BEGIN
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->destroy_datavector){
		free(this->datavector);
	}
	
	if(this->destroy_gammavector){
		free(this->gammavector);
	}

	if(this->destroy_thetavector){
		free(this->thetavector);
	}

	LOG_FUNC_END
}


/* print info about data */
void EdfData::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:           " << this->get_T() << std::endl;
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
