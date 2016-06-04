/** @file tsdata_global.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSDATAGLOBAL_H
#define	PASC_TSDATAGLOBAL_H

#ifndef USE_PETSCVECTOR
 #error 'TSDATA_GLOBAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "model/tsmodel_global.h"

namespace pascinference {

/* TSData GLOBAL for petsc */ 
class TSData_Global: public GeneralData {
	protected:
		TSModel_Global *tsmodel; /**< pointer to used time-series model on the data */

		GeneralVector<PetscVector> *datavector; /**< global vector with data of dimension based on model */
		bool destroy_datavector; /**< destroy datavector in destructor? if I am an owner, then TRUE */ 

		GeneralVector<PetscVector> *gammavector; /**< the characteristic functions of clustered models */
		bool destroy_gammavector;

		GeneralVector<PetscVector> *thetavector; /**< parameters of models */
		bool destroy_thetavector;

	public:
		TSData_Global(TSModel_Global &tsmodel);
		TSData_Global(TSModel_Global &tsmodel, GeneralVector<PetscVector> &datavector, GeneralVector<PetscVector> &gammavector, GeneralVector<PetscVector> &thetavector);
		~TSData_Global();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		void printcontent(ConsoleOutput &output) const;
		void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		std::string get_name() const;

		void cut_gamma() const;

		int get_T() const;
		int get_xdim() const;
		int get_K() const;

		/* GET functions */
		TSModel_Global *get_model() const;
		GeneralVector<PetscVector> *get_datavector() const;
		GeneralVector<PetscVector> *get_gammavector() const;
		GeneralVector<PetscVector> *get_thetavector() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* datavector is given */ //TODO: this is wrong
TSData_Global::TSData_Global(TSModel_Global &tsmodel, GeneralVector<PetscVector> &datavector, GeneralVector<PetscVector> &gammavector, GeneralVector<PetscVector> &thetavector){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* set new datavector */
	// TODO: control compatibility with this->tsmodel->get_datavectorlength();
	this->datavector = &datavector;
	destroy_datavector = false; /* this datavector is not my */

	/* set new gammavector and thetavector */
	// TODO: control compatibility with this->tsmodel->get_gammavectorlength(), this->tsmodel->get_thetavectorlength();
	this->gammavector = &gammavector;
	destroy_gammavector = false;

	this->thetavector = &thetavector;
	destroy_thetavector = false;

	LOG_FUNC_END
}

/* no datavector provided - prepare own data vector */
TSData_Global::TSData_Global(TSModel_Global &tsmodel){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new vectors based on model */
	Vec datavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
	TRY( VecSetSizes(datavector_Vec,this->tsmodel->get_datavectorlength_local(),this->tsmodel->get_datavectorlength_global()) );
	TRY( VecSetFromOptions(datavector_Vec) );
	this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
	this->destroy_datavector = true;

	Vec gammavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&gammavector_Vec) );
	TRY( VecSetSizes(gammavector_Vec,this->tsmodel->get_gammavectorlength_local(),this->tsmodel->get_gammavectorlength_global()) );
	TRY( VecSetFromOptions(gammavector_Vec) );
	this->gammavector = new GeneralVector<PetscVector>(gammavector_Vec);
	this->destroy_gammavector = true;

	Vec thetavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&thetavector_Vec) );
	TRY( VecSetSizes(thetavector_Vec,this->tsmodel->get_thetavectorlength_local(),this->tsmodel->get_thetavectorlength_global()) );
	TRY( VecSetFromOptions(thetavector_Vec) );
	this->thetavector = new GeneralVector<PetscVector>(thetavector_Vec);
	this->destroy_thetavector = true;

	LOG_FUNC_END
}

/* destructor */
TSData_Global::~TSData_Global(){
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
void TSData_Global::print(ConsoleOutput &output) const {
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
void TSData_Global::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
void TSData_Global::printcontent(ConsoleOutput &output) const {
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
void TSData_Global::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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

std::string TSData_Global::get_name() const {
	return "Time-series Data";
}



/* ---------- GET functions --------- */
int TSData_Global::get_T() const{
	if(this->tsmodel){
		return this->tsmodel->get_T();
	} else {
		return 0;
	}
}

int TSData_Global::get_xdim() const{
	if(this->tsmodel){
		return this->tsmodel->get_xdim();
	} else {
		return 0;
	}
}

int TSData_Global::get_K() const{
	if(this->tsmodel){
		return this->tsmodel->get_K();
	} else {
		return 0;
	}
}

TSModel_Global *TSData_Global::get_model() const{
	return this->tsmodel;
}

GeneralVector<PetscVector> *TSData_Global::get_datavector() const{
	return this->datavector;
}

GeneralVector<PetscVector> *TSData_Global::get_gammavector() const{
	return this->gammavector;
}

GeneralVector<PetscVector> *TSData_Global::get_thetavector() const{
	return this->thetavector;
}

void TSData_Global::cut_gamma() const{
	LOG_FUNC_BEGIN

	int max_id;
	double max_value;
	
	int K = get_K();
	int gamma_t = this->tsmodel->get_gammavectorlength_local()/(double)K;
	
	double *gamma_arr;
	TRY( VecGetArray(gammavector->get_vector(),&gamma_arr) );
	
	int t,k;
	for(t = 0; t < gamma_t; t++){
		/* find max value */
		max_id = 0;
		max_value = gamma_arr[t];
		for(k = 1; k < K; k++){
			if(gamma_arr[k*gamma_t + t] > max_value){
				max_id = k;
				max_value = gamma_arr[k*gamma_t + t];
			}
		}
		
		/* set new values */
		for(k = 0; k < K; k++){
			if(k == max_id){
				gamma_arr[k*gamma_t + t] = 1.0;
			} else {
				gamma_arr[k*gamma_t + t] = 0.0;
			}
		}


	}

	TRY( VecRestoreArray(gammavector->get_vector(),&gamma_arr) );
	
	LOG_FUNC_END
}


} /* end namespace */

#endif
