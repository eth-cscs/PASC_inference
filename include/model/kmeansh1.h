#ifndef KMEANSH1MODEL_H
#define	KMEANSH1MODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "model/tsmodel.h"
#include "algebra.h"

namespace pascinference {

/* KMEANSH1MODEL */ 
template<class VectorBase>
class KmeansH1Model: public TSModel<VectorBase> {
	protected:
		/* variables */
//		int T; /* length of time-series */
//		int K; /* number of clusters */
//		int dim; /* number of components in each time-step */

	public:
		KmeansH1Model(int T, int dim, int K);
		~KmeansH1Model();

		void print(std::ostream &output) const;

		int get_datavectorlength();
		int get_gammavectorlength();
		int get_thetavectorlength();
		
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
KmeansH1Model<VectorBase>::KmeansH1Model(int newT, int newdim, int newK) {
	if(DEBUG_MODE >= 100) std::cout << "(KmeansH1Model)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;

	/* prepare auxiliary vectors */

}

/* destructor */
template<class VectorBase>
KmeansH1Model<VectorBase>::~KmeansH1Model(){
	if(DEBUG_MODE >= 100) std::cout << "(KmeansH1Model)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about problem */
template<class VectorBase>
void KmeansH1Model<VectorBase>::print(std::ostream &output) const {
	output << " KmeansH1Model" << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:    " << this->T << std::endl;
	output << "  - K:    " << this->K << std::endl;
	output << "  - dim:  " << this->dim << std::endl;
		
}

/* get the appropriate length of datavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_datavectorlength(){
	return this->dim*this->T;
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_gammavectorlength(){
	return this->K*this->T;
}

/* get the appropriate length of thetavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_thetavectorlength(){
	return this->K*this->dim;
}


} /* end namespace */

#endif
