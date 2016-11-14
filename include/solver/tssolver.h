/** @file tssolver.h
 *  @brief Time-series solver.
 *
 *  @author Lukas Pospisil
 */
 
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

#define TSSOLVER_DEFAULT_MAXIT 300
#define TSSOLVER_DEFAULT_EPS 1e-9
#define TSSOLVER_DEFAULT_INIT_PERMUTE true

#define TSSOLVER_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class TSSolver
 *  \brief for solving time-series problems
 *
*/
template<class VectorBase>
class TSSolver: public GeneralSolver {
	protected:
		TSData<VectorBase> *tsdata; /**< tsdata on which the solver operates */

		GeneralSolver *gammasolver; /**< to solve inner gamma problem */
		bool gammasolved;			/**< when gamma is provided, this is true */

		GeneralSolver *thetasolver; /**< to solve inner theta problem */
		bool thetasolved;			/**< when theta is provided, this is true */
		
		TSModel<VectorBase> *model; /**< pointer to used time-series model */

		int it_sum; /**< sum of all iterations */
		int it_last; /**< number of interations in last solution */
		int annealing; /**< nuber of annealing steps */

		double L; /**< function value */
		double deltaL; /**< value of stopping criteria */
		std::ostringstream gammasolver_status; /**< status of gammasolver in the best annealing state */
		std::ostringstream thetasolver_status; /**< status of thetasolver in the best annealing state */
		std::ostringstream gammasolver_shortinfo_header; /**< shortinfo of gammasolver in the best annealing state */
		std::ostringstream gammasolver_shortinfo_values; /**< shortinfo of gammasolver in the best annealing state */
		std::ostringstream thetasolver_shortinfo_header; /**< shortinfo of thetasolver in the best annealing state */
		std::ostringstream thetasolver_shortinfo_values; /**< shortinfo of thetasolver in the best annealing state */

		Timer timer_solve; /**< total solution time of algorithm */
		Timer timer_gamma_solve; /**< timer for solving gamma problem */
		Timer timer_theta_solve; /**< timer for solving theta problem */
		Timer timer_gamma_update; /**< timer for updating gamma problem */
		Timer timer_theta_update; /**< timer for updating theta problem */

		bool init_permute;					/**< permute initial approximation or not */
		int debugmode;						/**< basic debug mode schema [0/1/2/3] */
		bool debug_print_annealing;			/**< print info about annealing steps */
		bool debug_print_it;				/**< print simple info about outer iterations */
		bool debug_print_theta;				/**< print theta solver info */
		bool debug_print_theta_solution;	/**< print solution of theta problem in each iteration */ 
		bool debug_print_gamma;				/**< print gamma solver info */
		bool debug_print_gamma_solution;	/**< print solution of gamma problem in each iteration */

		/* temp vectors for annealing */
		GeneralVector<VectorBase> *gammavector_temp;
		GeneralVector<VectorBase> *thetavector_temp;

		void prepare_temp_annealing();
		void destroy_temp_annealing();
		void set_settings_from_console();
		
		void gammavector_permute() const;
	public:
		TSSolver();
		TSSolver(TSData<VectorBase> &new_tsdata, int annealing=1);
		~TSSolver();

		virtual void solve();
		
		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		virtual void printstatus(ConsoleOutput &output) const;
		virtual void printstatus(std::ostringstream &output) const;
		virtual void printtimer(ConsoleOutput &output) const;
		virtual void printshort(std::ostringstream &header, std::ostringstream &values) const;
		virtual void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		virtual std::string get_name() const;

		virtual TSData<VectorBase> *get_data() const;

		virtual GeneralSolver *get_gammasolver() const;
		virtual GeneralSolver *get_thetasolver() const;

		void set_solution_theta(double *Theta);
		void set_annealing(int annealing);
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* set predefined debug mode schema */
template<class VectorBase>
void TSSolver<VectorBase>::set_settings_from_console(){
	LOG_FUNC_BEGIN
	
	consoleArg.set_option_value("tssolver_maxit", &this->maxit, TSSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("tssolver_eps", &this->eps, TSSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("tssolver_init_permute", &this->init_permute, TSSOLVER_DEFAULT_INIT_PERMUTE);

	/* set debug mode */
	consoleArg.set_option_value("tssolver_debugmode", &debugmode, TSSOLVER_DEFAULT_DEBUGMODE);

	debug_print_annealing = false;
	debug_print_it = false;
	debug_print_theta = false;
	debug_print_theta_solution = false; 
	debug_print_gamma = false;
	debug_print_gamma_solution = false; 

	if(debugmode == 1){
		debug_print_annealing = true;
	}

	if(debugmode == 2){
		debug_print_annealing = true;
		debug_print_it = true;
	}

	if(debugmode == 3){
		debug_print_annealing = true;
		debug_print_it = true;
		debug_print_theta = true;
		debug_print_gamma = true;
	}

	consoleArg.set_option_value("tssolver_debug_print_annealing", 		&debug_print_annealing, 		debug_print_annealing);
	consoleArg.set_option_value("tssolver_debug_print_it",				&debug_print_it, 				debug_print_it);
	consoleArg.set_option_value("tssolver_debug_print_theta", 			&debug_print_theta, 			debug_print_theta);
	consoleArg.set_option_value("tssolver_debug_print_theta_solution", 	&debug_print_theta_solution, 	debug_print_theta_solution);
	consoleArg.set_option_value("tssolver_debug_print_gamma", 			&debug_print_gamma, 			debug_print_gamma);
	consoleArg.set_option_value("tssolver_debug_print_gamma_solution", 	&debug_print_gamma_solution, 	debug_print_gamma_solution);

	LOG_FUNC_END
}

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
	this->annealing = 1;

	set_settings_from_console();

	this->L = std::numeric_limits<double>::max();
	this->deltaL = std::numeric_limits<double>::max();
	gammasolver_status.str("");
	thetasolver_status.str("");
	gammasolver_shortinfo_header.str("");
	gammasolver_shortinfo_values.str("");
	thetasolver_shortinfo_header.str("");
	thetasolver_shortinfo_values.str("");

	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();

	this->gammasolved = false;
	this->thetasolved = false;

	LOG_FUNC_END
}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(TSData<VectorBase> &new_tsdata, int annealing){
	LOG_FUNC_BEGIN

	/* set provided data and model */
	tsdata = &new_tsdata;
	model = tsdata->get_model(); 

	/* maybe the implementation of model is wrong, so I rather set it at first to NULL */
	gammasolver = NULL;
	thetasolver = NULL;

	/* we can initialize solvers - based on model */
	model->initialize_gammasolver(&gammasolver);
	model->initialize_thetasolver(&thetasolver);	

	/* set settings */
	set_settings_from_console();

	/* iteration counters */
	this->it_sum = 0;
	this->it_last = 0;
	this->annealing = annealing;

	/* initial value of object function */
	this->L = std::numeric_limits<double>::max();
	this->deltaL = std::numeric_limits<double>::max();
	gammasolver_status.str("");
	thetasolver_status.str("");
	gammasolver_shortinfo_header.str("");
	gammasolver_shortinfo_values.str("");
	thetasolver_shortinfo_header.str("");
	thetasolver_shortinfo_values.str("");

	/* initialize timers */
	this->timer_solve.restart();	
	this->timer_gamma_solve.restart();
	this->timer_theta_solve.restart();
	this->timer_gamma_update.restart();
	this->timer_theta_update.restart();

	this->gammasolved = false;
	this->thetasolved = false;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	LOG_FUNC_BEGIN

	/* destroy child solvers - based on model */
	if(gammasolver){
		model->finalize_gammasolver(&gammasolver);	
	}
	if(thetasolver){
		model->finalize_thetasolver(&thetasolver);	
	}

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* print settings */
	output <<  " - maxit:        " << this->maxit << std::endl;
	output <<  " - annealing:    " << this->annealing << std::endl;
	output <<  " - eps:          " << this->eps << std::endl;
	output <<  " - debugmode:   " << this->debugmode << std::endl;
	output <<  " - init_permute: " << this->init_permute << std::endl;

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

	output_global << this->get_name() << std::endl;

	/* print settings */
	output_global <<  " - maxit:        " << this->maxit << std::endl;
	output_global <<  " - eps:          " << this->eps << std::endl;
	output_global <<  " - debugmode:   " << this->debugmode << std::endl;
	output_global <<  " - init_permute: " << this->init_permute << std::endl;
	output_global <<  " - annealing:    " << this->annealing << std::endl;

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
	output <<  " - deltaL:      " << this->deltaL << std::endl;
	output <<  " - gammasolver_status:" << std::endl;
	output.push();
	output << gammasolver_status.str() << std::endl;
	output.pop();
	output <<  " - thetasolver_status:" << std::endl;
	output.push();
	output << thetasolver_status.str() << std::endl;
	output.pop();
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
	
	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - it:      " << std::setw(25) << this->it_last << std::endl;
	output <<  "      - L:       " << std::setw(25) << this->L << std::endl;
	output <<  "      - deltaL:  " << std::setw(25) << this->deltaL << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it all =          " << this->it_sum << std::endl;
	output <<  " - AIC =             " << tsdata->get_aic() << std::endl;
	output <<  " - annealing =       " << this->annealing << std::endl;
	output <<  " - init_permute =    " << this->init_permute << std::endl;
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
void TSSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	printshort_sum(header,values);
	
	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	header << "it all, ";
	values << this->it_sum << ", ";

	header << "t all, ";
	values << this->timer_solve.get_value_sum() << ", ";

	header << "t gamma update, ";
	values << this->timer_gamma_update.get_value_sum() << ", ";

	header << "t gamma solve, ";
	values << this->timer_gamma_solve.get_value_sum() << ", ";

	header << "t theta update, ";
	values << this->timer_theta_update.get_value_sum() << ", ";

	header << "t theta solve, ";
	values << this->timer_theta_solve.get_value_sum() << ", ";
	
	/* from best annealing step: */
	header << gammasolver_shortinfo_header.str();
	values << gammasolver_shortinfo_values.str();
	header << thetasolver_shortinfo_header.str();
	values << thetasolver_shortinfo_values.str();
	
	gammasolver->printshort_sum(header, values);
	thetasolver->printshort_sum(header, values);

	LOG_FUNC_END
}

template<class VectorBase>
std::string TSSolver<VectorBase>::get_name() const {
	return "TSSolver<VectorBase>";
}

template<class VectorBase>
void TSSolver<VectorBase>::set_annealing(int annealing) {
	this->annealing = annealing;
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

	/* now the gammasolver and thetasolver should be specified and prepared */

	/* variables */
	double L; /* object function value */
	double L_old;
	double deltaL;
	double aic;

	int it, it_annealing, it_gammasolver, it_thetasolver;

	/* prepare temp vectors for gamma and theta if there is more annealing steps */
	if(annealing > 1){
		prepare_temp_annealing();
	}

	/* annealing cycle */
	coutMaster.push();
	for(it_annealing=0;it_annealing < this->annealing;it_annealing++){
		if(debug_print_annealing){
			coutMaster <<  "- annealing = " << it_annealing << std::endl;
		}
		
		/* permute initial approximation subject to decomposition */
		if(this->init_permute){
			gammavector_permute();
		}
		
		/* couter for iterations inside annealing */
		it_gammasolver = 0;
		it_thetasolver = 0;

		/* initialize value of object function */
		L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
		deltaL = L;

		/* project init approximation to feasible set */
//		gammasolver->get_data()->get_feasibleset()->project(*gammadata->get_x0());
		
		/* main cycle */
		coutMaster.push();
		for(it=0;it < this->maxit;it++){
			if(debug_print_it){
				coutMaster <<  "it = " << it << std::endl;
			}

			/* --- COMPUTE Theta --- */
			if(!thetasolved){
				this->timer_theta_update.start();
				 model->update_thetasolver(thetasolver);
				this->timer_theta_update.stop();

				this->timer_theta_solve.start();
				 thetasolver->solve();
				this->timer_theta_solve.stop();
			}

			/* print info about theta solver */
			if(debug_print_theta){
				/* print info about cost function */
				coutMaster << " theta solver:" << std::endl;
				coutMaster << "  - ";
				coutMaster << "it = " << std::setw(6) << thetasolver->get_it() << ", ";
				coutMaster << "time_update = " << std::setw(12) << this->timer_theta_update.get_value_last() << ", ";
				coutMaster << "time_solve = " << std::setw(12) << this->timer_theta_solve.get_value_last() << std::endl;
				if(this->debugmode >= 10){
					coutMaster.push();
					thetasolver->printstatus(coutAll);
					coutMaster.pop();
				}
				coutAll.synchronize();

			}

			if(debug_print_theta_solution){
				coutMaster <<  "- thetasolver content:" << std::endl;
				coutMaster.push();
				 thetasolver->printcontent(coutMaster);
				coutMaster.pop();
			}

			/* --- COMPUTE gamma --- */
			if(!gammasolved){
				this->timer_gamma_update.start();
				 model->update_gammasolver(gammasolver);
				this->timer_gamma_update.stop();

				this->timer_gamma_solve.start();
				 gammasolver->solve();
				this->timer_gamma_solve.stop();
			}

			/* print info about gammasolver */
			if(debug_print_gamma){
				/* print info about cost function */
				coutMaster << " gamma solver:" << std::endl;
				coutMaster << "  - ";
				coutMaster << "it = " << std::setw(6) << gammasolver->get_it() << ", ";
				coutMaster << "time_update = " << std::setw(12) << this->timer_gamma_update.get_value_last() << ", ";
				coutMaster << "time_solve = " << std::setw(12) << this->timer_gamma_solve.get_value_last() << std::endl;
				if(this->debugmode >= 10){
					coutMaster.push();
					gammasolver->printstatus(coutAll);
					coutMaster.pop();
				}
				coutAll.synchronize();
			}
			if(debug_print_gamma_solution){
				coutMaster <<  "- gammasolver content:" << std::endl;
				coutMaster.push();
				gammasolver->printcontent(coutMaster);
				coutMaster.pop();
			}

			/* compute stopping criteria if the problem was not solved yet */
			L_old = L;
			L = model->get_L(gammasolver,thetasolver);
			deltaL = std::abs(L - L_old);

			/* log L */
			LOG_FX2(L,"L")
			LOG_FX2(deltaL,"deltaL")

			if(debug_print_it){
				/* print info about cost function */
				coutMaster << " outer loop status:" << std::endl;			
				coutMaster << "  - ";
				coutMaster << "L_old = " << std::setw(12) << L_old << ", ";
				coutMaster << "L = " << std::setw(12) << L << ", ";
				coutMaster << "|L - L_old| = " << std::setw(12) << deltaL << std::endl;
			}

			/* global stopping criteria */
			if(deltaL < this->eps){
				break;
			}

			/* update counter for outer annealing iterations */
			it_gammasolver += gammasolver->get_it();
			it_thetasolver += thetasolver->get_it();

		}
		coutMaster.pop();

		this->it_sum += it;
		this->it_last = it;

		/* compute AIC */ 
		aic = model->get_aic(L);

		if(debug_print_annealing){
			coutMaster << "  AIC=" << std::setw(12) << aic;
			coutMaster << ", L=" << std::setw(7) << L;
			coutMaster << ", it=" << std::setw(6) << it;
			coutMaster << ", it_gamma=" << std::setw(6) << it_gammasolver;
//		coutMaster << ", it_theta=" << std::setw(6) << it_thetasolver;
			coutMaster << std::endl;
		}

		/* if there is no other annealing steps, we are not using temp storage and store results directly */
		if((annealing <= 1) || (aic < tsdata->get_aic() && annealing > 1)){
			/* if this value is smaller then previous, then store it */
			if(aic < tsdata->get_aic() && annealing > 1){
				*gammavector_temp = *(tsdata->get_gammavector());
				*thetavector_temp = *(tsdata->get_thetavector());
			}

			tsdata->set_aic(aic);
			this->L = L;
			this->deltaL = deltaL;

			/* update status strings */
			gammasolver_status.str("");
			gammasolver->printstatus(gammasolver_status);
			thetasolver_status.str("");
			thetasolver->printstatus(thetasolver_status);
			/* update shortinfo strings */
			gammasolver_shortinfo_header.str("");
			gammasolver_shortinfo_values.str("");
			gammasolver->printshort(gammasolver_shortinfo_header,gammasolver_shortinfo_values);
			thetasolver_shortinfo_header.str("");
			thetasolver_shortinfo_values.str("");
			thetasolver->printshort(thetasolver_shortinfo_header,thetasolver_shortinfo_values);
		}

		/* if there are more annealing steps, then prepare new initial guess */
		if(it_annealing < this->annealing-1){
				tsdata->get_gammavector()->set_random();
		}
	}
	coutMaster.pop();
		
	/* destroy temp vectors for gamma and theta if there is more annealing steps */
	if(annealing > 1){
		*(tsdata->get_gammavector()) = *gammavector_temp;
		*(tsdata->get_thetavector()) = *thetavector_temp;

		destroy_temp_annealing();
	}

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

template<class VectorBase>
void TSSolver<VectorBase>::prepare_temp_annealing(){
	LOG_FUNC_BEGIN

	gammavector_temp = new GeneralVector<VectorBase>();
	thetavector_temp = new GeneralVector<VectorBase>();

	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::destroy_temp_annealing(){
	LOG_FUNC_BEGIN

	free(gammavector_temp);
	free(thetavector_temp);

	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::gammavector_permute() const{
	LOG_FUNC_BEGIN

	/* create new Vec and store into it permuted values */
	Vec gammavector_new;
	TRYCXX( VecDuplicate(tsdata->get_gammavector()->get_vector(), &gammavector_new) );
	this->tsdata->get_decomposition()->permute_TRK(tsdata->get_gammavector()->get_vector(), gammavector_new, false);

	/* copy values to original vector */
	TRYCXX( VecCopy(gammavector_new, tsdata->get_gammavector()->get_vector()) );

	/* destroy new vector */
	TRYCXX( VecDestroy(&gammavector_new) );

	LOG_FUNC_END
}

template<class VectorBase>
void TSSolver<VectorBase>::set_solution_theta(double *Theta) {
	LOG_FUNC_BEGIN

	double *theta_arr;
	TRYCXX(VecGetArray(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	for(int k=0;k<this->model->get_thetavectorlength_local();k++){
		theta_arr[k] = Theta[k];
	}
	TRYCXX(VecRestoreArray(tsdata->get_thetavector()->get_vector(), &theta_arr) );

	this->thetasolved = true;
	LOG_FUNC_END
}


}
} /* end namespace */

#endif
