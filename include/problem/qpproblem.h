#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalproblem.h"
#include "data/qpdata.h"
#include "solver/qpsolver.h"


namespace pascinference {

/* QPProblem */
template<class VectorBase>
class QPProblem: public GeneralProblem {
		QPData<VectorBase> *data; /* data of QP problem */
		QPSolver<VectorBase> *solver; /* solver of QP problem */
	
	public:
		QPProblem();
		~QPProblem();

		void print(std::ostream &output) const;
		void solve(); /* solve QPProblem */

		/* set data */
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
QPProblem<VectorBase>::QPProblem(){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)CONSTRUCTOR" << std::endl;
	
	/* create new data */
	this->data = new QPData<VectorBase>();
}

/* destructor */
template<class VectorBase>
QPProblem<VectorBase>::~QPProblem(){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)DESTRUCTOR" << std::endl;
	
	/* destroy data */
	free(this->data);
}


/* print info about problem */
template<class VectorBase>
void QPProblem<VectorBase>::print(std::ostream &output) const {
	output << "QPProblem" << std::endl;
	output << *data;
}

/* solve QPProblem */
template<class VectorBase>
void QPProblem<VectorBase>::solve(){
	std::cout << "We are solving the problem" << std::endl;
}

/* ----------------- SET ----------------- */

template<class VectorBase>
void QPProblem<VectorBase>::set_A(const GeneralMatrix<VectorBase> &newA){
	this->data->set_A(newA);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_b(const GeneralVector<VectorBase> &newb){
	this->data->set_b(newb);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_x(const GeneralVector<VectorBase> &newx){
	this->data->set_x(newx);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_x0(const GeneralVector<VectorBase> &newx0){
	this->data->set_x0(newx0);
}


} /* end namespace */

#endif
