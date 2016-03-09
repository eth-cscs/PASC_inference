#ifndef QPDATA_H
#define	QPDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "generalfeasibleset.h"
#include "algebra.h"

namespace pascinference {

/* QPData */ 
template<class VectorBase>
class QPData: public GeneralData {
	public:
		QPData();
		~QPData();

		void print(std::ostream &output) const;
		std::string get_name() const;

		/* variables */
		GeneralMatrix<VectorBase> *A; /* Hessian matrix */
		GeneralVector<VectorBase> *b; /* RHS vector, linear term */
		GeneralVector<VectorBase> *x0; /* initial approximation */
		
		GeneralFeasibleSet<VectorBase> *feasibleset; /* feasible set */

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
QPData<VectorBase>::QPData(){
	if(DEBUG_MODE >= 100) std::cout << "(QPData)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->A = NULL;
	this->b = NULL;
	this->x0 = NULL;
	this->feasibleset = NULL;

}

/* destructor */
template<class VectorBase>
QPData<VectorBase>::~QPData(){
	if(DEBUG_MODE >= 100) std::cout << "(QPData)DESTRUCTOR" << std::endl;
	
}


/* print info about problem */
template<class VectorBase>
void QPData<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - A:            ";
	if(this->A){
		output << "YES" << std::endl; // TODO: get matrix name 
	} else {
		output << "NO" << std::endl;
	}
	output << "  - b:            ";
	if(this->b){
		output << "YES (size: " << this->b->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - x0:           ";
	if(this->x0){
		output << "YES (size: " << this->x0->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - feasible_set: ";
	if(this->feasibleset){
		output << "YES (" << this->feasibleset->get_name() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
		
}

template<class VectorBase>
std::string QPData<VectorBase>::get_name() const {
	return "QP Data";
}

} /* end namespace */

#endif
