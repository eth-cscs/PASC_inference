#ifndef PASC_DIAGSOLVER_H
#define	PASC_DIAGSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/diagdata.h"

#define DIAGSOLVER_DEFAULT_DEBUG_MODE 0;

namespace pascinference {


/* DiagSolver */ 
template<class VectorBase>
class DiagSolver: public GeneralSolver {
	protected:
		Timer timer_solve; /* total solution time of SPG algorithm */
		Timer timer_dot; /* the sum of time necessary to compute dot_products */

		DiagData<VectorBase> *diagdata; /* data on which the solver operates */

	public:

		DiagSolver();
		DiagSolver(DiagData<VectorBase> &new_diagdata); 
		~DiagSolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void printstatus(ConsoleOutput &output) const;		
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		DiagData<VectorBase> *get_data() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(){
	LOG_FUNC_BEGIN
	
	diagdata = NULL;

	LOG_FUNC_END
}

template<class VectorBase>
DiagSolver<VectorBase>::DiagSolver(DiagData<VectorBase> &new_diagdata){
	LOG_FUNC_BEGIN

	diagdata = &new_diagdata;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
DiagSolver<VectorBase>::~DiagSolver(){
	LOG_FUNC_BEGIN
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void DiagSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debug_mode: " << this->debug_mode << std::endl;

	/* print data */
	if(diagdata){
		coutMaster.push();
		diagdata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}


template<class VectorBase>
void DiagSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void DiagSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =  " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_dot =    " << this->timer_dot.get_value_sum() << std::endl;

	LOG_FUNC_END
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
void DiagSolver<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); 

	this->timer_dot.start(); 
	 TRY( VecPointwiseDivide(diagdata->get_x()->get_vector(),diagdata->get_b()->get_vector(),diagdata->get_a()->get_vector() ) );
	this->timer_dot.stop(); 

	diagdata->get_x()->valuesUpdate();
	
	this->timer_solve.stop(); 

	LOG_FUNC_END
}

#endif

/* --------------- MINLIN ----------------- */
#ifdef USE_MINLIN

typedef minlin::threx::HostVector<double> MinlinHostVector;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;

template<>
void DiagSolver<MinlinHostVector>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); 

	typedef GeneralVector<MinlinHostVector> (&pVector);

	/* pointers to data */
	pVector a = *(diagdata->get_a());
	pVector b = *(diagdata->get_b());
	pVector x = *(diagdata->get_x());
	
	this->timer_dot.start(); 
	 int i;
	 for(i=0;i<x.size();i++){
		x(i) = b(i)/a(i);
	 }
	this->timer_dot.stop(); 
	
	this->timer_solve.stop(); 

	this->it = 1;
	this->fx = 0;

	/* write info to log file */
	LOG_IT(this->it)
	LOG_FX(this->fx)	

	LOG_FUNC_END
}

template<>
void DiagSolver<MinlinDeviceVector>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); 

	typedef GeneralVector<MinlinDeviceVector> (&pVector);

	/* pointers to data */
	pVector a = *(diagdata->get_a());
	pVector b = *(diagdata->get_b());
	pVector x = *(diagdata->get_x());
	
	this->timer_dot.start(); 
	 int i;
	 for(i=0;i<x.size();i++){
		x(i) = b(i)/a(i);
	 }
	this->timer_dot.stop(); 

	this->timer_solve.stop(); 

	LOG_FUNC_END
}

#endif

} /* end namespace */

#endif
