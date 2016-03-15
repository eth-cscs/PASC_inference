#ifndef QPRESULT_H
#define	QPRESULT_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalresult.h"
#include "algebra.h"

namespace pascinference {

/* QPResult */ 
template<class VectorBase>
class QPResult: public GeneralResult {
	public:
		QPResult();
		~QPResult();

		void print(std::ostream &output) const;
		std::string get_name() const;

		/* variables */
		GeneralVector<VectorBase> *x; /* solution */


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
QPResult<VectorBase>::QPResult(){
	if(DEBUG_MODE >= 100) coutMaster << "(QPResult)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->x = NULL;

}

/* destructor */
template<class VectorBase>
QPResult<VectorBase>::~QPResult(){
	if(DEBUG_MODE >= 100) coutMaster << "(QPResult)DESTRUCTOR" << std::endl;
	
}


/* print info about problem */
template<class VectorBase>
void QPResult<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - x:     ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
		
}

template<class VectorBase>
std::string QPResult<VectorBase>::get_name() const {
	return "QP Result";
}



} /* end namespace */

#endif
