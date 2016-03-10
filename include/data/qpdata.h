#ifndef PASC_QPDATA_H
#define	PASC_QPDATA_H

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
	private:
		/* variables */
		GeneralMatrix<VectorBase> *A; /* Hessian matrix */
		GeneralVector<VectorBase> *b; /* RHS vector, linear term */
		GeneralVector<VectorBase> *x0; /* initial approximation */
		GeneralVector<VectorBase> *x; /* solution */

		GeneralFeasibleSet<VectorBase> *feasibleset; /* feasible set */
	public:
		QPData();
		~QPData();

		void print(std::ostream &output) const;
		std::string get_name() const;

		/* set and get functions */
		void set_A(GeneralMatrix<VectorBase> *A) const;
		GeneralMatrix<VectorBase> *get_A() const;

		void set_b(GeneralVector<VectorBase> *b) const;
		GeneralVector<VectorBase> *get_b() const;

		void set_x0(GeneralVector<VectorBase> *x0) const;
		GeneralVector<VectorBase> *get_x0() const;

		void set_x(GeneralVector<VectorBase> *x) const;
		GeneralVector<VectorBase> *get_x() const;

		void set_feasibleset(GeneralFeasibleSet<VectorBase> *x0) const;
		GeneralFeasibleSet<VectorBase> *get_feasibleset() const;


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
	this->x = NULL;
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
	output << "  - x:           ";
	if(this->x0){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
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

/* ----- SET and GET functions --- */
template<class VectorBase>
void QPData<VectorBase>::set_A(GeneralMatrix<VectorBase> *A) const{
	this->A = A;
}

template<class VectorBase>
GeneralMatrix<VectorBase> *QPData<VectorBase>::get_A() const{
	return this->A;
}

template<class VectorBase>
void QPData<VectorBase>::set_b(GeneralVector<VectorBase> *b) const{
	this->b = b;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_b() const{
	return this->b;
}

template<class VectorBase>
void QPData<VectorBase>::set_x0(GeneralVector<VectorBase> *x0) const{
	this->x0 = x0;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_x0() const{
	return this->x0;
}

template<class VectorBase>
void QPData<VectorBase>::set_x(GeneralVector<VectorBase> *x) const{
	this->x = x;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_x() const{
	return this->x;
}

template<class VectorBase>
void QPData<VectorBase>::set_feasibleset(GeneralFeasibleSet<VectorBase> *feasibleset) const{
	this->feasibleset = feasibleset;
}

template<class VectorBase>
GeneralFeasibleSet<VectorBase> *QPData<VectorBase>::get_feasibleset() const{
	return this->feasibleset;
}


} /* end namespace */

#endif
