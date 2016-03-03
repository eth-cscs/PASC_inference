#ifndef TSMODEL_H
#define	TSMODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalmodel.h"
#include "algebra.h"

namespace pascinference {

/* Time-series MODEL */ 
template<class VectorBase>
class TSModel: public GeneralModel {
	protected:
		int T; /* length of time-series */
		int K; /* number of clusters */
		int dim; /* number of components in each time-step */

	public:
		TSModel();
		TSModel(int T, int dim, int K);
		~TSModel();

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		virtual int get_datavectorlength();
		virtual int get_gammavectorlength();
		virtual int get_thetavectorlength();

		int get_T() const;
		int get_K() const;
		int get_dim() const;
		
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSModel<VectorBase>::TSModel(){
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = 0;
	this->K = 0;
	this->dim = 0;
}


template<class VectorBase>
TSModel<VectorBase>::TSModel(int newT, int newdim, int newK){
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;
}

/* destructor */
template<class VectorBase>
TSModel<VectorBase>::~TSModel(){
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)DESTRUCTOR" << std::endl;
	
}


/* print info about model */
template<class VectorBase>
void TSModel<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:    " << this->T << std::endl;
	output << "  - K:    " << this->K << std::endl;
	output << "  - dim:  " << this->dim << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string TSModel<VectorBase>::get_name() const {
	return "General Time-Series Model";	
}

/* get the appropriate length of datavector */
template<class VectorBase>
int TSModel<VectorBase>::get_datavectorlength(){
	return 0;
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int TSModel<VectorBase>::get_gammavectorlength(){
	return 0;
}

/* get the appropriate length of thetavector */
template<class VectorBase>
int TSModel<VectorBase>::get_thetavectorlength(){
	return 0; 
}


/* ---- GET FUNCTIONS ---- */
template<class VectorBase>
int TSModel<VectorBase>::get_T() const{
	return T; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_dim() const{
	return dim; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_K() const{
	return K; 
}


} /* end namespace */

#endif
