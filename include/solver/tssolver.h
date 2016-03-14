#ifndef PASC_TSSOLVER_H
#define	PASC_TSSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/tsdata.h"
#include "model/tsmodel.h"


//temp
#include "solver/qpsolver.h"
#include "solver/diagsolver.h"


#define TSSOLVER_DEFAULT_MAXIT 1000;
#define TSSOLVER_DEFAULT_EPS 0.001;
#define TSSOLVER_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class TSSolverSetting : public GeneralSolverSetting {
	protected:

	public:
		TSSolverSetting() {
			this->maxit = TSSOLVER_DEFAULT_MAXIT;
			this->eps = TSSOLVER_DEFAULT_EPS;
			this->debug_mode = TSSOLVER_DEFAULT_DEBUG_MODE;
		};
		~TSSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << " TSSolverSettings:" << std::endl;
			output << "   - debug mode: " << this->debug_mode << std::endl;
			output << "   - maxit: " << this->maxit << std::endl;
			output << "   - eps: " << this->eps << std::endl;

		};
		
};


/* TSSolver */ 
template<class VectorBase>
class TSSolver: public GeneralSolver {
	protected:
		TSData<VectorBase> *tsdata; /* tsdata on which the solver operates */

		GeneralSolver *gammasolver; /* to solve inner gamma problem */
		GeneralSolver *thetasolver; /* to solve inner theta problem */

		TSModel<VectorBase> *model; /* pointer to used time-series model */

		int it; /**< actual iteration */

	public:
		TSSolverSetting setting;

		TSSolver();
		TSSolver(TSData<VectorBase> &new_tsdata); 
		//TODO: write constructor with given gammasolver and theta solver
		~TSSolver();

		/* with given gamma and theta solvertype */
		void solve(SolverType gammasolvertype, SolverType thetasolvertype);
		
		/* with one solvertype */
		virtual void solve(SolverType solvertype) {
			this->solve(solvertype,solvertype); 
		};
		
		/* without given solvertype */
		virtual void solve() {
			this->solve(SOLVER_AUTO); 
		};

		void print(std::ostream &output) const;
		void printstatus(std::ostream &output) const;
		void printtimer(std::ostream &output) const;
		std::string get_name() const;

		TSData<VectorBase> *get_data() const;
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSSolver<VectorBase>::TSSolver(){
	if(setting.debug_mode >= 100) std::cout << "(TSSolver)CONSTRUCTOR" << std::endl;
	
	tsdata = NULL;
	model = NULL;
	gammasolver = NULL; /* in this time, we don't know how to solve the problem */
	thetasolver = NULL; /* in this time, we don't know how to solve the problem */
	
}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(TSData<VectorBase> &new_tsdata){
	tsdata = &new_tsdata;
	model = tsdata->get_model(); 

	/* we can initialize solvers - based on model */
	model->initialize_gammasolver(&gammasolver, tsdata);	
	model->initialize_thetasolver(&thetasolver, tsdata);	
	
}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	if(setting.debug_mode >= 100) std::cout << "(TSSolver)DESTRUCTOR" << std::endl;

	/* destroy child solvers - based on model */
	if(gammasolver){
		model->finalize_gammasolver(&gammasolver, tsdata);	
	}
	if(thetasolver){
		model->finalize_thetasolver(&thetasolver, tsdata);	
	}
	
}


/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) std::cout << this->get_name() << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;
	
	/* if child solvers are specified, then print also info about it */	
	output << " Gamma Solver" << std::endl;
	if(gammasolver){
		output << *gammasolver << std::endl;
	} else {
		output << " - not set" << std::endl;
	}

	output << " Theta Solver" << std::endl;
	if(thetasolver){
		output << *thetasolver << std::endl;
	} else {
		output << " - not set" << std::endl;
	}


}

template<class VectorBase>
void TSSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)FUNCTION: printstatus" << std::endl;

	output << this->get_name() << std::endl;
	output << " - it:          " << this->it << std::endl;
	output << " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
}

template<class VectorBase>
void TSSolver<VectorBase>::printtimer(std::ostream &output) const {
	output << this->get_name() << std::endl;
	output << " - it =        " << this->it << std::endl;
	output << " - timers" << std::endl;

	output << " Gamma Solver" << std::endl;
	if(gammasolver){
		gammasolver->printtimer(output);
	} else {
		output << " - not set" << std::endl;
	}

	output << " Theta Solver" << std::endl;
	if(thetasolver){
		thetasolver->printtimer(output);
	} else {
		output << " - not set" << std::endl;
	}
	
}


template<class VectorBase>
std::string TSSolver<VectorBase>::get_name() const {
	return "Time-Series Solver";
}

/* solve the problem */
template<class VectorBase>
void TSSolver<VectorBase>::solve(SolverType gammasolvertype, SolverType thetasolvertype) {
	if(setting.debug_mode >= 100) std::cout << "(TSSolver)FUNCTION: solve" << std::endl;

	/* the gamma or theta solver wasn't specified yet */
	if(!gammasolver || !thetasolver){
		// TODO: give error - actually, gamma and theta solvers are created during constructor with tsdata - now we don't have tsdata, there is nothing to solve
	}

	/* which specific solver we can use to solve the problem? */
	if(gammasolvertype == SOLVER_AUTO){
			// TODO: here write sofisticated decision tree - maybe based on MODEL?
	} 
	if(thetasolvertype == SOLVER_AUTO){
			// TODO: here write sofisticated decision tree - maybe based on MODEL?
	} 

	/* now the gammasolver and thetasolver should be specified and prepared */

	/* variables */
	double L, L_old, deltaL; /* object function value */

	/* initialize value of object function */
	L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
	
	int it; 

	/* main cycle */
	for(it=0;it < setting.maxit;it++){
		Message_info_value(" - it = ",it);

		/* --- COMPUTE Theta --- */
		model->update_thetasolver(thetasolver, tsdata);
		thetasolver->solve(thetasolvertype);

		if(setting.debug_mode >= 10){
			thetasolver->printstatus(std::cout);
		}

		/* print Theta */
		std::cout << "Theta: " << *(tsdata->get_thetavector()) << std::endl;


		/* --- COMPUTE gamma --- */
		model->update_gammasolver(gammasolver, tsdata);
		gammasolver->solve(gammasolvertype);

		if(setting.debug_mode >= 10){
			gammasolver->printstatus(std::cout);
		}


		/* compute stopping criteria */
		L_old = L;
		L = model->get_L(gammasolver,thetasolver,tsdata);
		deltaL = std::abs(L - L_old);

		/* print info about cost function */
		Message_info_value("  - L_old       = ",L_old);
		Message_info_value("  - L           = ",L);
		Message_info_value("  - |L - L_old| = ",deltaL);

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < setting.eps){
			break;
		}
		
	}

	this->it = it;
	
}

template<class VectorBase>
TSData<VectorBase> *TSSolver<VectorBase>::get_data() const {
	return tsdata;
}



} /* end namespace */

#endif
