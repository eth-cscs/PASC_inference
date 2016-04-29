#ifndef PASC_TSSOLVER_GLOBAL_H
#define	PASC_TSSOLVER_GLOBAL_H

#ifndef USE_PETSCVECTOR
 #error 'TSSOLVER_GLOBAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/tsdata_global.h"
#include "model/tsmodel_global.h"


//temp
#include "solver/qpsolver.h"
#include "solver/diagsolver.h"


#define TSSOLVER_GLOBAL_DEFAULT_MAXIT 1000;
#define TSSOLVER_GLOBAL_DEFAULT_EPS 0.001;
#define TSSOLVER_GLOBAL_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class TSSolverGlobalSetting : public GeneralSolverSetting {
	protected:

	public:
		SolverType gammasolvertype;
		SolverType thetasolvertype;

		TSSolverGlobalSetting() {
			this->maxit = TSSOLVER_GLOBAL_DEFAULT_MAXIT;
			this->eps = TSSOLVER_GLOBAL_DEFAULT_EPS;
			this->debug_mode = TSSOLVER_GLOBAL_DEFAULT_DEBUG_MODE;

			this->gammasolvertype = SOLVER_AUTO;
			this->thetasolvertype = SOLVER_AUTO;
		};
		~TSSolverGlobalSetting() {};

		void print(std::ostream &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - debug mode: " << this->debug_mode << std::endl;
			output <<  " - maxit: " << this->maxit << std::endl;
			output <<  " - eps: " << this->eps << std::endl;

		};

		std::string get_name() const {
			return "Time-Serie-Global Solver Setting";
		};
		
		
};


/* TSSolver_Global_Global */ 
class TSSolver_Global: public GeneralSolver {
	protected:
		TSData_Global *tsdata; /**< tsdata on which the solver operates */

		GeneralSolver *gammasolver; /**< to solve inner gamma problem */
		GeneralSolver *thetasolver; /**< to solve inner theta problem */

		TSModel_Global *model; /**< pointer to used time-series model */

		int it_sum; /**< sum of all iterations */
		int it_last; /**< number of interations in last solution */

		double L; /**< function value */

		Timer timer_solve; /**< total solution time of algorithm */
		Timer timer_gamma_solve; /**< timer for solving gamma problem */
		Timer timer_theta_solve; /**< timer for solving theta problem */
		Timer timer_gamma_update; /**< timer for updating gamma problem */
		Timer timer_theta_update; /**< timer for updating theta problem */

	public:
		TSSolverGlobalSetting setting;

		TSSolver_Global();
		TSSolver_Global(TSData_Global &new_tsdata); 
		~TSSolver_Global();

		virtual void solve();
		
		void print(std::ostream &output) const;
		void print(std::ostream &output_global, std::ostream &output_local) const;
		
		void printstatus(std::ostream &output) const;
		void printtimer(std::ostream &output) const;
		std::string get_name() const;

		TSData_Global *get_data() const;
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
TSSolver_Global::TSSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(TSSolver_Global)CONSTRUCTOR" << std::endl;
	
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

TSSolver_Global::TSSolver_Global(TSData_Global &new_tsdata){
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
TSSolver_Global::~TSSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(TSSolver_Global)DESTRUCTOR" << std::endl;

	/* destroy child solvers - based on model */
	if(gammasolver){
		model->finalize_gammasolver(&gammasolver, tsdata);	
	}
	if(thetasolver){
		model->finalize_thetasolver(&thetasolver, tsdata);	
	}
	
}


/* print info about problem */
void TSSolver_Global::print(std::ostream &output) const {

	output <<  this->get_name() << std::endl;

	coutMaster.push();
	setting.print(output);
	coutMaster.pop();

	/* print data */
	if(tsdata){
		coutMaster.push();
		tsdata->print(output);
		coutMaster.pop();
	}

	/* if child solvers are specified, then print also info about it */	
	output <<  " Gamma Solver" << std::endl;
	if(gammasolver){
		coutMaster.push();
		gammasolver->print(output);
		coutMaster.pop();
	} else {
		output <<  " - not set" << std::endl;
	}

	output <<  " Theta Solver" << std::endl;
	if(thetasolver){
		coutMaster.push();
		thetasolver->print(output);
		coutMaster.pop();
	} else {
		output <<  " - not set" << std::endl;
	}
}

/* print info about problem */
void TSSolver_Global::print(std::ostream &output_global, std::ostream &output_local) const {

	output_global <<  this->get_name() << std::endl;

	coutMaster.push();
	setting.print(output_global);
	coutMaster.pop();

	/* print data */
	if(tsdata){
		coutMaster.push();
		tsdata->print(output_global, output_local);
		coutMaster.pop();
	}

	/* if child solvers are specified, then print also info about it */	
	output_global <<  " Gamma Solver" << std::endl;
	if(gammasolver){
		coutMaster.push();
		gammasolver->print(output_global);
		coutMaster.pop();
	} else {
		output_global <<  " - not set" << std::endl;
	}

	output_global <<  " Theta Solver" << std::endl;
	if(thetasolver){
		coutMaster.push();
		thetasolver->print(output_global);
		coutMaster.pop();
	} else {
		output_global <<  " - not set" << std::endl;
	}
}


void TSSolver_Global::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(SPGQPSolver)FUNCTION: printstatus" << std::endl;

	output <<  this->get_name() << std::endl;
	output <<  " - it:          " << this->it_last << std::endl;
	output <<  " - L:           " << this->L << std::endl;	
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
}

void TSSolver_Global::printtimer(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	output <<  " - it all =          " << this->it_sum << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =        " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_gamma_update = "  << this->timer_gamma_update.get_value_sum() << std::endl;
	output <<  "  - t_gamma_solve =  "  << this->timer_gamma_solve.get_value_sum() << std::endl;
	output <<  "  - t_theta_update = " << this->timer_theta_update.get_value_sum() << std::endl;
	output <<  "  - t_theta_solve =  " << this->timer_theta_solve.get_value_sum() << std::endl;

	output <<  " Gamma Solver" << std::endl;
	if(gammasolver){
		coutMaster.push();
		gammasolver->printtimer(output);
		coutMaster.pop();
	} else {
		output << " - not set" << std::endl;
	}

	output <<  " Theta Solver" << std::endl;
	if(thetasolver){
		coutMaster.push();
		thetasolver->printtimer(output);
		coutMaster.pop();
	} else {
		output << " - not set" << std::endl;
	}
	
}


std::string TSSolver_Global::get_name() const {
	return "Time-Series Solver";
}

/* solve the problem */
void TSSolver_Global::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(TSSolver_Global)FUNCTION: solve" << std::endl;

	this->timer_solve.start(); /* stop this timer in the end of solution */

	/* the gamma or theta solver wasn't specified yet */
	if(!gammasolver || !thetasolver){
		// TODO: give error - actually, gamma and theta solvers are created during constructor with tsdata - now we don't have tsdata, there is nothing to solve
		coutMaster << "Warning: gammasolver of thetasolver is not set yet!" << std::endl;
	}

	/* update settings of child solvers */ //TODO: this is not working at all
	gammasolver->setting.debug_mode = setting.debug_mode;
	thetasolver->setting.debug_mode = setting.debug_mode;

	/* now the gammasolver and thetasolver should be specified and prepared */

	/* variables */
	double L; /* object function value */
//	double L_old, deltaL;

	/* initialize value of object function */
	L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
	
	int it; 

	/* main cycle */
	coutMaster.push();
	for(it=0;it < setting.maxit;it++){
		coutMaster <<  "it = " << it << std::endl;

		/* --- COMPUTE Theta --- */
		coutAll << "------------------------ update theta solver" << std::endl;
		this->timer_theta_update.start();
		 model->update_thetasolver(thetasolver, tsdata);
		this->timer_theta_update.stop();

		/* barrier  - everything has to be synchronized now */
		TRY(PetscBarrier(NULL));

		coutAll << "------------------------ solve theta problem" << std::endl;
		this->timer_theta_solve.start();
//		 thetasolver->solve();
		this->timer_theta_solve.stop();

		/* print info about theta solver */
		if(setting.debug_mode >= 10){
			coutMaster <<  "- thetasolver status:" << std::endl;
			coutMaster.push();
			thetasolver->printstatus(coutMaster);
			coutMaster.pop();
		}
		if(setting.debug_mode >= 100){
			coutMaster <<  "- thetasolver info:" << std::endl;
			coutMaster.push();
			thetasolver->print(coutMaster);
			coutMaster.pop();
		}
		if(setting.debug_mode >= 101){
			coutMaster <<  "- thetasolver content:" << std::endl;
			coutMaster.push();
			thetasolver->printcontent(coutMaster);
			coutMaster.pop();
		}


		/* --- COMPUTE gamma --- */
		this->timer_gamma_update.start();
//		 model->update_gammasolver(gammasolver, tsdata);
		this->timer_gamma_update.stop();

		this->timer_gamma_solve.start();
//		 gammasolver->solve();
		this->timer_gamma_solve.stop();

		/* print info about gammasolver */
		if(setting.debug_mode >= 10){
			coutMaster <<  "- gammasolver status:" << std::endl;
			coutMaster.push();
			gammasolver->printstatus(coutMaster);
			coutMaster.pop();
		}
		if(setting.debug_mode >= 100){
			coutMaster <<  "- gammasolver info:" << std::endl;
			coutMaster.push();
			gammasolver->print(coutMaster);
			coutMaster.pop();
		}
		if(setting.debug_mode >= 101){
			coutMaster <<  "- gammasolver content:" << std::endl;
			coutMaster.push();
			gammasolver->printcontent(coutMaster);
			coutMaster.pop();
		}

		/* compute stopping criteria */
//		L_old = L;
//		L = model->get_L(gammasolver,thetasolver,tsdata);
//		deltaL = std::abs(L - L_old);

		/* print info about cost function */
//		coutMaster <<  " - L_old       = " << L_old << std::endl;
//		coutMaster <<  " - L           = " << L << std::endl;
//		coutMaster <<  " - |L - L_old| = " << deltaL << std::endl;

		/* end the main cycle if the change of function value is sufficient */
//		if (deltaL < setting.eps){
//			break;
//		}
		
	}
	coutMaster.pop();

	this->it_sum += it;
	this->it_last = it;
	this->L = L;

	this->timer_solve.stop(); /* stop this timer in the end of solution */
	
}

TSData_Global *TSSolver_Global::get_data() const {
	return tsdata;
}



} /* end namespace */

#endif
