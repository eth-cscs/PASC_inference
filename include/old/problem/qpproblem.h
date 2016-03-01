#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalproblem.h"
#include "data/qpdata.h"
#include "solver/qpsolver.h"
#include "result/qpresult.h"


namespace pascinference {

/* QPProblem */
template<class VectorBase>
class QPProblem: public GeneralProblem {
	protected:
		QPData<VectorBase> *data; /* data of QP problem */
		QPSolver<VectorBase> *solver; /* solver of QP problem */
		QPResult<VectorBase> *result; /* results of QP problem */
	
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
	
	/* create new data and results */
	this->data = new QPData<VectorBase>();
	this->result = new QPResult<VectorBase>();

	/* create new solver */
	this->solver = new QPSolver<VectorBase>(this->data,this->result);

}

/* destructor */
template<class VectorBase>
QPProblem<VectorBase>::~QPProblem(){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)DESTRUCTOR" << std::endl;
	
	/* destroy data */
	free(this->data);

	/* destroy result */
	free(this->result);

	/* destroy solver */
	free(this->solver);

}


/* print info about problem */
template<class VectorBase>
void QPProblem<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)OPERATOR <<" << std::endl;

	output << "QPProblem" << std::endl;
	output << *data;
	output << *solver;
	output << *result;
}

/* solve QPProblem */
template<class VectorBase>
void QPProblem<VectorBase>::solve(){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)FUNCTION solve" << std::endl;

	std::cout << "We are solving the problem" << std::endl;
}

/* ----------------- SET ----------------- */

template<class VectorBase>
void QPProblem<VectorBase>::set_A(const GeneralMatrix<VectorBase> &newA){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)FUNCTION set_A" << std::endl;

	this->data->set_A(newA);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_b(const GeneralVector<VectorBase> &newb){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)FUNCTION set_b" << std::endl;

	this->data->set_b(newb);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_x(const GeneralVector<VectorBase> &newx){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)FUNCTION set_x" << std::endl;

	this->result->set_x(newx);
}

template<class VectorBase>
void QPProblem<VectorBase>::set_x0(const GeneralVector<VectorBase> &newx0){
	if(DEBUG_MODE >= 100) std::cout << "(QPProblem)FUNCTION set_x0" << std::endl;

	this->data->set_x0(newx0);
}


} /* end namespace */

#endif
