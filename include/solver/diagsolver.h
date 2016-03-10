#ifndef PASC_DIAGSOLVER_H
#define	PASC_DIAGSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/diagdata.h"


namespace pascinference {

/* settings */
class DiagSolverSetting : public GeneralSolverSetting {
	protected:
		int maxit;
		double eps;

	public:
		DiagSolverSetting() {
		};
		~DiagSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  DiagSolverSettings:" << std::endl;
		};
		
};


/* DiagSolver */ 
template<class VectorBase>
class DiagSolver: public GeneralSolver {
	protected:
		const DiagData<VectorBase> *data; /* data on which the solver operates */

	public:
		DiagSolverSetting setting;

		DiagSolver();
		DiagSolver(const DiagData<VectorBase> &new_data); 
		~DiagSolver();

		virtual void solve(SolverType solvertype);
		virtual void solve() { this->solve(SOLVER_AUTO); };

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)CONSTRUCTOR" << std::endl;
	
	data = NULL;
}

template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(const DiagData<VectorBase> &new_data){
	data = &new_data;
}

/* destructor */
template<class VectorBase>
DiagSolver<VectorBase>::~DiagSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)DESTRUCTOR" << std::endl;
}


/* print info about problem */
template<class VectorBase>
void DiagSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)FUNCTION: print" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;

	/* print data */
	output << *data;
	
}

template<class VectorBase>
std::string DiagSolver<VectorBase>::get_name() const {
	return "Diagonal Solver";
}

/* solve the problem */
template<class VectorBase>
void DiagSolver<VectorBase>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)FUNCTION: solve" << std::endl;

	typedef GeneralVector<VectorBase> (&pVector);

	/* pointers to data */
	pVector a = *(data->get_a());
	pVector b = *(data->get_b());
	pVector x = *(data->get_x());
	
	x = a/b;
	
}


} /* end namespace */

#endif
