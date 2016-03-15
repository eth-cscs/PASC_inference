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

		void print(std::ostream &output) const {
			output << offset << this->get_name() << std::endl;
			output << offset << " - debug mode: " << this->debug_mode << std::endl;
			output << offset << " - maxit: " << this->maxit << std::endl;
			output << offset << " - eps: " << this->eps << std::endl;

		};

		std::string get_name() const {
			return "Time-Serie Solver Setting";
		};
		
		
};


/* TSSolver */ 
template<class VectorBase>
class TSSolver: public GeneralSolver {
	protected:
		TSData<VectorBase> *tsdata; /**< tsdata on which the solver operates */

		GeneralSolver *gammasolver; /**< to solve inner gamma problem */
		GeneralSolver *thetasolver; /**< to solve inner theta problem */

		TSModel<VectorBase> *model; /**< pointer to used time-series model */

		int it_sum; /**< sum of all iterations */
		int it_last; /**< number of interations in last solution */

		double L; /**< function value */

		Timer timer_solve; /**< total solution time of algorithm */
		Timer timer_gamma_solve; /**< timer for solving gamma problem */
		Timer timer_theta_solve; /**< timer for solving theta problem */
		Timer timer_gamma_update; /**< timer for updating gamma problem */
		Timer timer_theta_update; /**< timer for updating theta problem */

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
	if(setting.debug_mode >= 100) coutMaster << offset <<"(TSSolver)CONSTRUCTOR" << std::endl;
	
	tsdata = NULL;
	model = NULL;
	gammasolver = NULL; /* in this time, we don't know how to solve the problem */
	thetasolver = NULL; /* in this time, we don't know how to solve the problem */

	this->it_sum = 0;
	this->it_last = 0;

	this->L = std::numeric_limits<double>::max();

	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();

}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(TSData<VectorBase> &new_tsdata){
	tsdata = &new_tsdata;
	model = tsdata->get_model(); 

	/* we can initialize solvers - based on model */
	model->initialize_gammasolver(&gammasolver, tsdata);	
	model->initialize_thetasolver(&thetasolver, tsdata);	
	
	this->it_sum = 0;
	this->it_last = 0;

	this->L = std::numeric_limits<double>::max();

	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();
}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	if(setting.debug_mode >= 100) coutMaster << offset <<"(TSSolver)DESTRUCTOR" << std::endl;

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

	output << offset << this->get_name() << std::endl;

	offset.push();
	setting.print(output);
	offset.pop();

	/* print data */
	if(tsdata){
		offset.push();
		tsdata->print(output);
		offset.pop();
	}

	/* if child solvers are specified, then print also info about it */	
	output << offset << " Gamma Solver" << std::endl;
	if(gammasolver){
		offset.push();
		gammasolver->print(output);
		offset.pop();
	} else {
		output << offset << " - not set" << std::endl;
	}

	output << offset << " Theta Solver" << std::endl;
	if(thetasolver){
		offset.push();
		thetasolver->print(output);
		offset.pop();
	} else {
		output << offset << " - not set" << std::endl;
	}
}

template<class VectorBase>
void TSSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << offset <<"(SPGQPSolver)FUNCTION: printstatus" << std::endl;

	output << offset << this->get_name() << std::endl;
	output << offset << " - it:          " << this->it_last << std::endl;
	output << offset << " - L:           " << this->L << std::endl;	
	output << offset << " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
}

template<class VectorBase>
void TSSolver<VectorBase>::printtimer(std::ostream &output) const {
	output << offset << this->get_name() << std::endl;
	output << offset << " - it all =          " << this->it_sum << std::endl;
	output << offset << " - timers" << std::endl;
	output << offset << "  - t_solve =        " << this->timer_solve.get_value_sum() << std::endl;
	output << offset << "  - t_gamma_update = "  << this->timer_gamma_update.get_value_sum() << std::endl;
	output << offset << "  - t_gamma_solve =  "  << this->timer_gamma_solve.get_value_sum() << std::endl;
	output << offset << "  - t_theta_update = " << this->timer_theta_update.get_value_sum() << std::endl;
	output << offset << "  - t_theta_solve =  " << this->timer_theta_solve.get_value_sum() << std::endl;

	output << offset << " Gamma Solver" << std::endl;
	if(gammasolver){
		offset.push();
		gammasolver->printtimer(output);
		offset.pop();
	} else {
		output << " - not set" << std::endl;
	}

	output << offset << " Theta Solver" << std::endl;
	if(thetasolver){
		offset.push();
		thetasolver->printtimer(output);
		offset.pop();
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
	if(setting.debug_mode >= 100) coutMaster << offset <<"(TSSolver)FUNCTION: solve" << std::endl;

	this->timer_solve.start(); /* stop this timer in the end of solution */

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
	offset.push();
	for(it=0;it < setting.maxit;it++){
		coutMaster << offset << "it = " << it << std::endl;

		/* --- COMPUTE Theta --- */
		this->timer_theta_update.start();
		 model->update_thetasolver(thetasolver, tsdata);
		this->timer_theta_update.stop();

		this->timer_theta_solve.start();
		 thetasolver->solve(thetasolvertype);
		this->timer_theta_solve.stop();

		if(setting.debug_mode >= 10){
			offset.push();
			thetasolver->printstatus(coutMaster);
			offset.pop();
		}

		/* --- COMPUTE gamma --- */
		this->timer_gamma_update.start();
		 model->update_gammasolver(gammasolver, tsdata);
		this->timer_gamma_update.stop();

		this->timer_gamma_solve.start();
		 gammasolver->solve(gammasolvertype);
		this->timer_gamma_solve.stop();

		if(setting.debug_mode >= 10){
			offset.push();
			gammasolver->printstatus(coutMaster);
			offset.pop();
		}


		/* compute stopping criteria */
		L_old = L;
		L = model->get_L(gammasolver,thetasolver,tsdata);
		deltaL = std::abs(L - L_old);

		/* print info about cost function */
		coutMaster << offset << " - L_old       = " << L_old << std::endl;
		coutMaster << offset << " - L           = " << L << std::endl;
		coutMaster << offset << " - |L - L_old| = " << deltaL << std::endl;

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < setting.eps){
			break;
		}
		
	}
	offset.pop();

	this->it_sum += it;
	this->it_last = it;
	this->L = L;

	this->timer_solve.stop(); /* stop this timer in the end of solution */
	
}

template<class VectorBase>
TSData<VectorBase> *TSSolver<VectorBase>::get_data() const {
	return tsdata;
}



} /* end namespace */

#endif
