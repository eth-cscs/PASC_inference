#ifndef PASC_TSSOLVER_H
#define	PASC_TSSOLVER_H

#ifndef USE_PETSCVECTOR
 #error 'TSSOLVER is for PETSCVECTOR'
 typedef petscvector::PetscVector PetscVector;
#endif


#include "pascinference.h"

#include "data/tsdata.h"
#include "model/tsmodel.h"

//temp
#include "solver/qpsolver.h"
#include "solver/diagsolver.h"
#include "solver/qpsolver.h"

#define TSSOLVER_DEFAULT_MAXIT 1000
#define TSSOLVER_DEFAULT_EPS 0.001
#define TSSOLVER_DEFAULT_DEBUG_MODE 0

namespace pascinference {

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
		TSSolver();
		TSSolver(TSData<VectorBase> &new_tsdata); 
		~TSSolver();

		virtual void solve();
		
		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printstatus(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		TSData<VectorBase> *get_data() const;

		GeneralSolver *get_gammasolver() const;
		GeneralSolver *get_thetasolver() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSSolver<VectorBase>::TSSolver(){
	LOG_FUNC_BEGIN
	
	this->tsdata = NULL;
	this->model = NULL;
	this->gammasolver = NULL; /* in this time, we don't know how to solve the problem */
	this->thetasolver = NULL; /* in this time, we don't know how to solve the problem */

	this->it_sum = 0;
	this->it_last = 0;

	consoleArg.set_option_value("tssolver_maxit", &this->maxit, TSSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("tssolver_eps", &this->eps, TSSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("tssolver_debug_mode", &this->debug_mode, TSSOLVER_DEFAULT_DEBUG_MODE);

	this->L = std::numeric_limits<double>::max();

	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();

	LOG_FUNC_END
}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(TSData<VectorBase> &new_tsdata){
	LOG_FUNC_BEGIN

	/* set provided data and model */
	tsdata = &new_tsdata;
	model = tsdata->get_model(); 

	/* we can initialize solvers - based on model */
	model->initialize_gammasolver(&gammasolver, tsdata);	
	model->initialize_thetasolver(&thetasolver, tsdata);	

	/* set settings */
	consoleArg.set_option_value("tssolver_maxit", &this->maxit, TSSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("tssolver_eps", &this->eps, TSSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("tssolver_debug_mode", &this->debug_mode, TSSOLVER_DEFAULT_DEBUG_MODE);

	/* iteration counters */
	this->it_sum = 0;
	this->it_last = 0;

	/* initial value of object function */
	this->L = std::numeric_limits<double>::max();

	/* initialize timers */
	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	LOG_FUNC_BEGIN

	/* destroy child solvers - based on model */
	if(gammasolver){
		model->finalize_gammasolver(&gammasolver, tsdata);	
	}
	if(thetasolver){
		model->finalize_thetasolver(&thetasolver, tsdata);	
	}
	
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debug_mode: " << this->debug_mode << std::endl;

	/* print data */
	if(tsdata){
		output.push();
		tsdata->print(output);
		output.pop();
	}

	/* if child solvers are specified, then print also info about it */	
	output <<  " Gamma Solver" << std::endl;
	if(gammasolver){
		output.push();
		gammasolver->print(output);
		output.pop();
	} else {
		output <<  " - not set" << std::endl;
	}

	output <<  " Theta Solver" << std::endl;
	if(thetasolver){
		output.push();
		thetasolver->print(output);
		output.pop();
	} else {
		output <<  " - not set" << std::endl;
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << std::endl;
	output_global <<  " - eps:        " << this->eps << std::endl;
	output_global <<  " - debug_mode: " << this->debug_mode << std::endl;

	/* print data */
	if(tsdata){
		output_global.push();
		tsdata->print(output_global, output_local);
		output_global.pop();
	}

	/* if child solvers are specified, then print also info about it */	
	output_global <<  " Gamma Solver" << std::endl;
	if(gammasolver){
		output_global.push();
		gammasolver->print(output_global, output_local);
		output_global.pop();
	} else {
		output_global <<  " - not set" << std::endl;
	}

	output_global <<  " Theta Solver" << std::endl;
	if(thetasolver){
		output_global.push();
		thetasolver->print(output_global, output_local);
		output_global.pop();
	} else {
		output_global <<  " - not set" << std::endl;
	}
	
	output_local.synchronize();
	output_global.synchronize();

	LOG_FUNC_END
}


template<class VectorBase>
void TSSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it:          " << this->it_last << std::endl;
	output <<  " - L:           " << this->L << std::endl;	
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
	
	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

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

	LOG_FUNC_END
}


template<class VectorBase>
std::string TSSolver<VectorBase>::get_name() const {
	return "TSSolver<VectorBase>";
}

/* solve the problem */
template<class VectorBase>
void TSSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); /* stop this timer in the end of solution */

	std::ostringstream temp_ostream;
	
	/* the gamma or theta solver wasn't specified yet */
	if(!gammasolver || !thetasolver){
		// TODO: give error - actually, gamma and theta solvers are created during constructor with tsdata - now if we don't have tsdata, there is nothing to solve
		coutMaster << "Warning: gammasolver or thetasolver is not set yet!" << std::endl;
	}

	/* update settings of child solvers */ //TODO: this is not working at all
	gammasolver->debug_mode = this->debug_mode;

	/* now the gammasolver and thetasolver should be specified and prepared */

	/* variables */
	double L; /* object function value */
	double L_old;
	double deltaL;

	/* initialize value of object function */
	L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
	deltaL = L;
	
	int it; 

	/* main cycle */
	coutMaster.push();
	for(it=0;it < this->maxit;it++){
		coutMaster <<  "it = " << it << std::endl;

		/* --- COMPUTE Theta --- */
		this->timer_theta_update.start();
		 model->update_thetasolver(thetasolver, tsdata);
		this->timer_theta_update.stop();

		/* barrier  - everything has to be synchronized now */
		TRY(PetscBarrier(NULL));

		this->timer_theta_solve.start();
		 thetasolver->solve();
		this->timer_theta_solve.stop();

		TRY(PetscBarrier(NULL));

		/* print info about theta solver */
		if(this->debug_mode >= 2){
			/* print info about cost function */
			coutMaster << " theta solver:" << std::endl;
			coutAll << "  - ";
			coutAll << "it = " << std::setw(6) << thetasolver->get_it() << ", ";
			coutAll << "time_update = " << std::setw(12) << this->timer_theta_update.get_value_last() << ", ";
			coutAll << "time_solve = " << std::setw(12) << this->timer_theta_solve.get_value_last() << std::endl;
			if(this->debug_mode >= 10){
				coutMaster.push();
				thetasolver->printstatus(coutAll);
				coutMaster.pop();
			}
			coutAll.synchronize();

		}
					
		if(this->debug_mode >= 100){
			coutMaster <<  "- thetasolver info:" << std::endl;
			coutMaster.push();
			thetasolver->print(coutMaster);
			coutMaster.pop();
		}
		if(this->debug_mode >= 101){
			coutMaster <<  "- thetasolver content:" << std::endl;
			coutMaster.push();
			thetasolver->printcontent(coutMaster);
			coutMaster.pop();
		}


		/* --- COMPUTE gamma --- */
		this->timer_gamma_update.start();
		 model->update_gammasolver(gammasolver, tsdata);
		this->timer_gamma_update.stop();

		TRY(PetscBarrier(NULL));

		this->timer_gamma_solve.start();
		 gammasolver->solve();
		this->timer_gamma_solve.stop();

		/* print info about gammasolver */
		if(this->debug_mode >= 2){
			/* print info about cost function */
			coutMaster << " gamma solver:" << std::endl;
			coutAll << "  - ";
			coutAll << "it = " << std::setw(6) << gammasolver->get_it() << ", ";
			coutAll << "time_update = " << std::setw(12) << this->timer_gamma_update.get_value_last() << ", ";
			coutAll << "time_solve = " << std::setw(12) << this->timer_gamma_solve.get_value_last() << std::endl;
			if(this->debug_mode >= 10){
				coutMaster.push();
				gammasolver->printstatus(coutAll);
				coutMaster.pop();
			}
			coutAll.synchronize();
		}
		if(this->debug_mode >= 100){
			coutMaster <<  "- gammasolver info:" << std::endl;
			coutMaster.push();
			gammasolver->print(coutMaster);
			coutMaster.pop();
		}
		if(this->debug_mode >= 101){
			coutMaster <<  "- gammasolver content:" << std::endl;
			coutMaster.push();
			gammasolver->printcontent(coutMaster);
			coutMaster.pop();
		}

		/* compute stopping criteria if the problem was not solved yet */
		L_old = L;
		L = model->get_L(gammasolver,thetasolver,tsdata);
		deltaL = std::abs(L - L_old);

		/* log L */
		LOG_FX2(L,"L")
		LOG_FX2(deltaL,"deltaL")

		if(this->debug_mode >= 2){
			/* print info about cost function */
			coutMaster << " outer loop status:" << std::endl;			
			coutAll << "  - ";
			coutAll << "L_old = " << std::setw(12) << L_old << ", ";
			coutAll << "L = " << std::setw(12) << L << ", ";
			coutAll << "|L - L_old| = " << std::setw(12) << deltaL << std::endl;
			coutAll.synchronize();
		}

		/* global stopping criteria */
		if(deltaL < this->eps){
			break;
		}

	}
	coutMaster.pop();

	this->it_sum += it;
	this->it_last = it;
	this->L = L;

	this->timer_solve.stop(); /* stop this timer in the end of solution */

	LOG_IT(this->it_last)
	LOG_FX(this->L)
	
	LOG_FUNC_END
}

template<class VectorBase>
TSData<VectorBase> *TSSolver<VectorBase>::get_data() const {
	return tsdata;
}

template<class VectorBase>
GeneralSolver *TSSolver<VectorBase>::get_gammasolver() const {
	return gammasolver;
}

template<class VectorBase>
GeneralSolver *TSSolver<VectorBase>::get_thetasolver() const {
	return thetasolver;
}


} /* end namespace */

#endif