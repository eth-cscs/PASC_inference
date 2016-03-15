#ifndef TSRESULT_H
#define	TSRESULT_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalresult.h"
#include "model/tsmodel.h"
#include "algebra.h"

namespace pascinference {

/* maybe TSModel is not defined yet. */
template<class VectorBase> class TSModel;

/* QPResult */ 
template<class VectorBase>
class TSResult: public GeneralResult {
	protected:
		TSModel<VectorBase> *tsmodel; /* pointer to used time-series model on the data */

		GeneralVector<VectorBase> *gammavector; /* the characteristic functions of clustered models */
		GeneralVector<VectorBase> *thetavector; /* the characteristic functions of clustered models */
		bool destroy_vectors; /* destroy gamma & theta vectors in the destructor? */

	public:
		TSResult();
		TSResult(TSModel<VectorBase> &tsmodel);
		TSResult(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &gammavector, GeneralVector<VectorBase> &thetavector);
		~TSResult();

		void print(std::ostream &output) const;
		std::string get_name() const;

		int get_T() const;
		int get_dim() const;
		int get_K() const;

		TSModel<VectorBase> *get_model() const;

		GeneralVector<VectorBase> *get_gammavector() const;
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSResult<VectorBase>::TSResult(){
	if(DEBUG_MODE >= 100) coutMaster << "(QPResult)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->tsmodel = NULL;

	this->gammavector = NULL;
	this->thetavector = NULL;
	destroy_vectors = false;

}

/* gammavector and thetavector are given */
template<class VectorBase>
TSResult<VectorBase>::TSResult(TSModel<VectorBase> &tsmodel, GeneralVector<VectorBase> &gammavector, GeneralVector<VectorBase> &thetavector){
	if(DEBUG_MODE >= 100) coutMaster << "(TSResult)CONSTRUCTOR model, gammavector, thetavector" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* set new datavector */
	// TODO: control compatibility with this->tsmodel->get_gammavectorlength(), this->tsmodel->get_thetavectorlength();
	this->gammavector = &gammavector;
	this->thetavector = &thetavector;
	destroy_vectors = false; /* these vectors are not my */
		
}

/* no datavector provided - prepare own data vector */
template<class VectorBase>
TSResult<VectorBase>::TSResult(TSModel<VectorBase> &tsmodel){
	if(DEBUG_MODE >= 100) coutMaster << "(TSResult)CONSTRUCTOR model" << std::endl;

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new result vectors */
	int gammavector_length = this->tsmodel->get_gammavectorlength();
	this->gammavector = new GeneralVector<VectorBase>(gammavector_length);

	int thetavector_length = this->tsmodel->get_thetavectorlength();
	this->thetavector = new GeneralVector<VectorBase>(thetavector_length);

	destroy_vectors = true;
	
}


/* destructor */
template<class VectorBase>
TSResult<VectorBase>::~TSResult(){
	if(DEBUG_MODE >= 100) coutMaster << "(QPResult)DESTRUCTOR" << std::endl;
	
	/* if I created a result vectors, then I should also be able to destroy them */
	if(this->destroy_vectors){
		free(this->gammavector);
		free(this->thetavector);
	}
	
}


/* print info about problem */
template<class VectorBase>
void TSResult<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:           " << this->get_T() << std::endl;
	output << "  - dim:         " << this->get_dim() << std::endl;
	output << "  - K:           " << this->get_K() << std::endl;
	output << "  - model:       " << this->tsmodel->get_name() << std::endl;
	output << "  - gammavector: ";
	if(this->gammavector){
		output << "YES (size: " << this->gammavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - thetavector: ";
	if(this->thetavector){
		output << "YES (size: " << this->thetavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
}

template<class VectorBase>
std::string TSResult<VectorBase>::get_name() const {
	return "Time-series Result";
}



/* ---------- GET functions --------- */
template<class VectorBase>
int TSResult<VectorBase>::get_T() const{
	if(this->tsmodel){
		return this->tsmodel->get_T();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSResult<VectorBase>::get_dim() const{
	if(this->tsmodel){
		return this->tsmodel->get_dim();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSResult<VectorBase>::get_K() const{
	if(this->tsmodel){
		return this->tsmodel->get_K();
	} else {
		return 0;
	}
}

template<class VectorBase>
TSModel<VectorBase> *TSResult<VectorBase>::get_model() const{
	return this->tsmodel;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSResult<VectorBase>::get_gammavector() const{
	return this->gammavector;
}


} /* end namespace */

#endif
