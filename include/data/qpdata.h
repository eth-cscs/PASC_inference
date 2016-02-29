#ifndef QPDATA_H
#define	QPDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "algebra.h"

namespace pascinference {

/* QPData */ 
template<class VectorBase>
class QPData: public GeneralData {
	private:
		const GeneralMatrix<VectorBase> *A; /* Hessian matrix */
		const GeneralVector<VectorBase> *b; /* RHS vector, linear term */
		const GeneralVector<VectorBase> *x0; /* initial approximation */
		const GeneralVector<VectorBase> *x; /* solution */
		
	public:
		QPData();
		~QPData();

		void print(std::ostream &output) const;
		
		void set_A(const GeneralMatrix<VectorBase> &newA); 
		void set_b(const GeneralVector<VectorBase> &newb); 
		void set_x0(const GeneralVector<VectorBase> &newx0); 
		void set_x(const GeneralVector<VectorBase> &newx); 

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

}

/* destructor */
template<class VectorBase>
QPData<VectorBase>::~QPData(){
	if(DEBUG_MODE >= 100) std::cout << "(QPData)DESTRUCTOR" << std::endl;
	
}


/* print info about problem */
template<class VectorBase>
void QPData<VectorBase>::print(std::ostream &output) const {
	output << " QPData" << std::endl;
	
	/* give information about presence of the data */
	output << "  - A:     ";
	if(this->A){
		output << "YES" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - b:     ";
	if(this->b){
		output << "YES" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - x0:    ";
	if(this->x0){
		output << "YES" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - x:     ";
	if(this->x){
		output << "YES" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
		
}

/* ------------- SET ----------------- */
template<class VectorBase>
void QPData<VectorBase>::set_A(const GeneralMatrix<VectorBase> &newA){
	this->A = &newA;
}

template<class VectorBase>
void QPData<VectorBase>::set_b(const GeneralVector<VectorBase> &newb){
	this->b = &newb;
}

template<class VectorBase>
void QPData<VectorBase>::set_x(const GeneralVector<VectorBase> &newx){
	this->x = &newx;
}

template<class VectorBase>
void QPData<VectorBase>::set_x0(const GeneralVector<VectorBase> &newx0){
	this->x0 = &newx0;
}


} /* end namespace */

#endif
