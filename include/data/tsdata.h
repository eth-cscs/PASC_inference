#ifndef TSDATA_H
#define	TSDATA_H

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
		TSModel<VectorBase> *tsmodel; /* pointer to used time-series model on the data */

		GeneralVector<VectorBase> *datavector; /* vector with data of dimension based on model */
		bool destroy_datavector; /* destroy datavector in destructor? if I am an owner, then TRUE */ 

	public:
		TSData();
		TSData(TSModel<VectorBase> &tsmodel);
		TSData(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &datavector);
		~TSData();

		void print(std::ostream &output) const;
		std::string get_name() const;

		int get_T() const;
		int get_dim() const;
		
		TSModel<VectorBase> *get_model() const;
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSData<VectorBase>::TSData(){
	if(DEBUG_MODE >= 100) std::cout << "(TSData)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->model = NULL;

	this->datavector = NULL;
	destroy_datavector = false;
}

/* datavector is given */
template<class VectorBase>
TSData<VectorBase>::TSData(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &datavector){
	if(DEBUG_MODE >= 100) std::cout << "(TSData)CONSTRUCTOR model, datavector" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* set new datavector */
	// TODO: control compatibility with this->tsmodel->get_datavectorlength();
	this->datavector = &datavector;
	destroy_datavector = false; /* this datavector is not my */
		
}

/* no datavector provided - prepare own data vector */
template<class VectorBase>
TSData<VectorBase>::TSData(TSModel<VectorBase> &tsmodel){
	if(DEBUG_MODE >= 100) std::cout << "(TSData)CONSTRUCTOR model" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new datavector */
	int datavector_length = this->tsmodel->get_datavectorlength();
	this->datavector = new GeneralVector<VectorBase>(datavector_length);
	destroy_datavector = true;
}

/* destructor */
template<class VectorBase>
TSData<VectorBase>::~TSData(){
	if(DEBUG_MODE >= 100) std::cout << "(TSData)DESTRUCTOR" << std::endl;
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->destroy_datavector){
		free(this->datavector);
	}
}


/* print info about problem */
template<class VectorBase>
void TSData<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:          " << this->get_T() << std::endl;
	output << "  - dim:        " << this->get_dim() << std::endl;
	output << "  - model:      " << this->tsmodel->get_name() << std::endl;
	output << "  - datavector: ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
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
TSModel<VectorBase> *TSData<VectorBase>::get_model() const{
	return this->tsmodel;
}

} /* end namespace */

#endif
