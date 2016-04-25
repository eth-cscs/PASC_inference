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
#include "generaldata.h"
#include "algebra.h"
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

		GeneralVector<PetscVector> *u; /**< external influence */
		bool destroy_u;

	public:
		TSData_Global();
		TSData_Global(TSModel_Global &tsmodel);
		TSData_Global(TSModel_Global &tsmodel, GeneralVector<PetscVector> &datavector, GeneralVector<PetscVector> &gammavector, GeneralVector<PetscVector> &thetavector, GeneralVector<PetscVector> &u);
		~TSData_Global();

		void print(std::ostream &output) const;
		void printcontent(std::ostream &output_global, std::ostream &output_local) const;
		std::string get_name() const;

		int get_T() const;
		int get_xdim() const;
		int* get_K() const;

		/* GET functions */
		TSModel_Global *get_model() const;
		GeneralVector<PetscVector> *get_datavector() const;
		GeneralVector<PetscVector> *get_u() const;
		GeneralVector<PetscVector> *get_gammavector() const;
		GeneralVector<PetscVector> *get_thetavector() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
TSData_Global::TSData_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->tsmodel = NULL;

	this->datavector = NULL;
	destroy_datavector = false;

	this->u = NULL;
	destroy_u = false;

	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

}

/* datavector is given */ //TODO: this is wrong
TSData_Global::TSData_Global(TSModel_Global &tsmodel, GeneralVector<PetscVector> &datavector, GeneralVector<PetscVector> &gammavector, GeneralVector<PetscVector> &thetavector, GeneralVector<PetscVector> &u){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData)CONSTRUCTOR model, datavector, gammavector, thetavector, u" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* set new datavector */
	// TODO: control compatibility with this->tsmodel->get_datavectorlength();
	this->datavector = &datavector;
	destroy_datavector = false; /* this datavector is not my */

	this->u = NULL;
	destroy_u = false;

	/* set new gammavector and thetavector */
	// TODO: control compatibility with this->tsmodel->get_gammavectorlength(), this->tsmodel->get_thetavectorlength();
	this->gammavector = &gammavector;
	destroy_gammavector = false;

	this->thetavector = &thetavector;
	destroy_thetavector = false;

}

/* no datavector provided - prepare own data vector */
TSData_Global::TSData_Global(TSModel_Global &tsmodel){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData_Global)CONSTRUCTOR model" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new vectors based on model */
	Vec datavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
	TRY( VecSetSizes(datavector_Vec,this->tsmodel->get_datavectorlength_local(),this->tsmodel->get_datavectorlength_global()) );
	TRY( VecSetFromOptions(datavector_Vec) );
	this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
	this->destroy_datavector = true;

	if(this->tsmodel->get_ulength_global() > 0){
		Vec u_Vec;
		TRY( VecCreate(PETSC_COMM_WORLD,&u_Vec) );
		TRY( VecSetSizes(u_Vec,this->tsmodel->get_ulength_local(),this->tsmodel->get_ulength_global()) );
		TRY( VecSetFromOptions(u_Vec) );
		this->u = new GeneralVector<PetscVector>(u_Vec);
		this->destroy_u = true;
	} else {
		this->u = NULL;
		this->destroy_u = false;
	}

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

}

/* destructor */
TSData_Global::~TSData_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData_Global)DESTRUCTOR" << std::endl;
	
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

	if(this->destroy_u){
		free(this->u);
	}

}


/* print info about data */
void TSData_Global::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:           " << this->get_T() << std::endl;
	output <<  " - xdim:        " << this->get_xdim() << std::endl;
	output <<  " - K:           ";
	print_array(output, this->tsmodel->get_Knum(), this->tsmodel->get_K());
	output << std::endl;
	output <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	output <<  " - datavector:  ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<  " - u:           ";
	if(this->u){
		output << "YES (size: " << this->u->size() << ")" << std::endl;
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
	
}

/* print content of all data */
void TSData_Global::printcontent(std::ostream &output_global,std::ostream &output_local) const {
	output_global <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output_local <<  " - datavector: ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local <<  " - u:          ";
	if(this->u){
		output_local << *this->u << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local <<  " - gammavector: ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local <<  " - thetavector: ";
	if(this->thetavector){
		output_local << *this->thetavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
		
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

int* TSData_Global::get_K() const{
	if(this->tsmodel){
		return this->tsmodel->get_K();
	} else {
		return NULL;
	}
}

TSModel_Global *TSData_Global::get_model() const{
	return this->tsmodel;
}

GeneralVector<PetscVector> *TSData_Global::get_datavector() const{
	return this->datavector;
}

GeneralVector<PetscVector> *TSData_Global::get_u() const{
	return this->u;
}

GeneralVector<PetscVector> *TSData_Global::get_gammavector() const{
	return this->gammavector;
}

GeneralVector<PetscVector> *TSData_Global::get_thetavector() const{
	return this->thetavector;
}


} /* end namespace */

#endif
