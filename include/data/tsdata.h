#ifndef PASC_TSDATA_H
#define	PASC_TSDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "algebra.h"
#include "model/tsmodel.h"


namespace pascinference {

/* TSData */ 
template<class VectorBase>
class TSData: public GeneralData {
	protected:
		TSModel<VectorBase> *tsmodel; /**< pointer to used time-series model on the data */

		GeneralVector<VectorBase> *datavector; /**< vector with data of dimension based on model */
		bool destroy_datavector; /**< destroy datavector in destructor? if I am an owner, then TRUE */ 

		GeneralVector<VectorBase> *gammavector; /**< the characteristic functions of clustered models */
		bool destroy_thetavector;

		GeneralVector<VectorBase> *thetavector; /**< parameters of models */
		bool destroy_gammavector;

		GeneralVector<VectorBase> *u; /**< external influence */
		bool destroy_u;

	public:
		TSData();
		TSData(TSModel<VectorBase> &tsmodel);
		TSData(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &datavector, GeneralVector<VectorBase> &gammavector, GeneralVector<VectorBase> &thetavector, GeneralVector<VectorBase> &u);
		~TSData();

		void print(std::ostream &output) const;
		void printcontent(std::ostream &output) const;
		std::string get_name() const;

		int get_T() const;
		int get_dim() const;
		int get_K() const;

		/* GET functions */
		TSModel<VectorBase> *get_model() const;
		GeneralVector<VectorBase> *get_datavector() const;
		GeneralVector<VectorBase> *get_u() const;
		GeneralVector<VectorBase> *get_gammavector() const;
		GeneralVector<VectorBase> *get_thetavector() const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSData<VectorBase>::TSData(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->model = NULL;

	this->datavector = NULL;
	destroy_datavector = false;

	this->u = NULL;
	destroy_u = false;

	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

}

/* datavector is given */
template<class VectorBase>
TSData<VectorBase>::TSData(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &datavector, GeneralVector<VectorBase> &gammavector, GeneralVector<VectorBase> &thetavector, GeneralVector<VectorBase> &u){
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
template<class VectorBase>
TSData<VectorBase>::TSData(TSModel<VectorBase> &tsmodel){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData)CONSTRUCTOR model" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new datavectors */
	int datavector_length = this->tsmodel->get_datavectorlength();
	this->datavector = new GeneralVector<VectorBase>(datavector_length);
	destroy_datavector = true;

	int u_length = this->tsmodel->get_ulength();
	if(u_length > 0){
		this->u = new GeneralVector<VectorBase>(u_length);
		destroy_u = true;
	} else {
		this->u = NULL;
		destroy_u = false;
	}
	
	int gammavector_length = this->tsmodel->get_gammavectorlength();
	this->gammavector = new GeneralVector<VectorBase>(gammavector_length);
	destroy_gammavector = true;

	int thetavector_length = this->tsmodel->get_thetavectorlength();
	this->thetavector = new GeneralVector<VectorBase>(thetavector_length);
	destroy_thetavector = true;
	
	
}

/* destructor */
template<class VectorBase>
TSData<VectorBase>::~TSData(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSData)DESTRUCTOR" << std::endl;
	
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
template<class VectorBase>
void TSData<VectorBase>::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:           " << this->get_T() << std::endl;
	output <<  " - dim:         " << this->get_dim() << std::endl;
	output <<  " - K:           " << this->get_K() << std::endl;
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
template<class VectorBase>
void TSData<VectorBase>::printcontent(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output <<  " - datavector: ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - u:          ";
	if(this->u){
		output << *this->u << std::endl;
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

		
}

template<class VectorBase>
std::string TSData<VectorBase>::get_name() const {
	return "Time-series Data";
}



/* ---------- GET functions --------- */
template<class VectorBase>
int TSData<VectorBase>::get_T() const{
	if(this->tsmodel){
		return this->tsmodel->get_T();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_dim() const{
	if(this->tsmodel){
		return this->tsmodel->get_dim();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_K() const{
	if(this->tsmodel){
		return this->tsmodel->get_K();
	} else {
		return 0;
	}
}

template<class VectorBase>
TSModel<VectorBase> *TSData<VectorBase>::get_model() const{
	return this->tsmodel;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_datavector() const{
	return this->datavector;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_u() const{
	return this->u;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_gammavector() const{
	return this->gammavector;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_thetavector() const{
	return this->thetavector;
}


} /* end namespace */

#endif
