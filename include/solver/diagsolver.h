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
		DiagData<VectorBase> *diagdata; /* data on which the solver operates */

	public:
		DiagSolverSetting setting;

		DiagSolver();
		DiagSolver(DiagData<VectorBase> &new_diagdata); 
		~DiagSolver();

		virtual void solve(SolverType solvertype);
		virtual void solve() { this->solve(SOLVER_AUTO); };

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		virtual DiagData<VectorBase> *get_data() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)CONSTRUCTOR" << std::endl;
	
	diagdata = NULL;
}

template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(DiagData<VectorBase> &new_diagdata){
	diagdata = &new_diagdata;
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
	output << *diagdata;
	
}

template<class VectorBase>
std::string DiagSolver<VectorBase>::get_name() const {
	return "Diagonal Solver";
}

template<class VectorBase>
DiagData<VectorBase> *DiagSolver<VectorBase>::get_data() const {
	return diagdata;
}

/* ---------- PETSCVECTOR ------------ */
#ifdef USE_PETSCVECTOR
typedef petscvector::PetscVector PetscVector;

/* Petsc: constructor from given right PetscVector */
template<>
void DiagSolver<PetscVector>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)FUNCTION: solve" << std::endl;

	TRY( VecPointwiseDivide(diagdata->get_x()->get_vector(),diagdata->get_b()->get_vector(),diagdata->get_a()->get_vector() ) );

	diagdata->get_x()->valuesUpdate();
	
}

#endif

/* --------------- MINLIN ----------------- */
#ifdef USE_MINLIN

typedef minlin::threx::HostVector<double> MinlinHostVector;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;

template<>
void DiagSolver<MinlinHostVector>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)FUNCTION: solve" << std::endl;

	typedef GeneralVector<MinlinHostVector> (&pVector);

	/* pointers to data */
	pVector a = *(diagdata->get_a());
	pVector b = *(diagdata->get_b());
	pVector x = *(diagdata->get_x());
	
	int i;
	for(i=0;i<x.size();i++){
		x(i) = b(i)/a(i);
	}
	
}

template<>
void DiagSolver<MinlinDeviceVector>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(DiagSolver)FUNCTION: solve" << std::endl;

	typedef GeneralVector<MinlinDeviceVector> (&pVector);

	/* pointers to data */
	pVector a = *(diagdata->get_a());
	pVector b = *(diagdata->get_b());
	pVector x = *(diagdata->get_x());
	
	int i;
	for(i=0;i<x.size();i++){
		x(i) = b(i)/a(i);
	}
	
}

#endif

} /* end namespace */

#endif
