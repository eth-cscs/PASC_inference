/** @file entropysolvernewton.h
 *  @brief Solver which solves problem with integrals using Newton method
 *  
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERNEWTON_H
#define	PASC_ENTROPYSOLVERNEWTON_H

#include "general/solver/generalsolver.h"
#include "general/data/entropydata.h"

/* include integration algorithms */
#include "general/algebra/integration/entropyintegrationdlib.h"
#include "general/algebra/integration/entropyintegrationcuba.h"

#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT 1000
#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT_AXB 100
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS 1e-6
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS_AXB 1e-6
#define ENTROPYSOLVERNEWTON_DEFAULT_INTEGRATION_EPS 1e-10
#define ENTROPYSOLVERNEWTON_DEFAULT_INTEGRATION_TYPE 0
#define ENTROPYSOLVERNEWTON_DEFAULT_NEWTON_COEFF 0.9
#define ENTROPYSOLVERNEWTON_DEFAULT_DEBUGMODE 0

#define ENTROPYSOLVERNEWTON_MONITOR false

namespace pascinference {
namespace solver {

/** \class EntropySolverNewton
 *  \brief Solver which solves problem with integrals using SPG algorithm
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolverNewton: public GeneralSolver {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		double *fxs;							/**< function values in clusters */
		double *gnorms;							/**< norm of gradient in clusters (stopping criteria) */
	
		int maxit_Axb;							/**< maximum number of iterations for solving inner system */
		double eps_Axb;							/**< the precision of KSP solver */
		int *it_sums;							/**< sums of all iterations for each cluster */
		int *it_lasts;							/**< number of iterations from last solve() call for each cluster */
		int *itAxb_sums;						/**< sums of all cg iterations for each cluster */
		int *itAxb_lasts;						/**< sums of all cg iterations in this outer iteration */
		double newton_coeff;					/**< newton step-size coefficient x_{k+1} = x_k + coeff*delta */

		EntropyIntegration<VectorBase> *entropyintegration;	/**< instance of integration tool */
		
		Timer timer_compute_moments;			/**< time for computing moments from data */
		Timer timer_solve; 						/**< total solution time of Newton algorithm */
		Timer timer_Axb; 						/**< total solution time of KSP algorithm */
		Timer timer_update; 					/**< total time of vector updates */
		Timer timer_g;			 				/**< total time of gradient computation */
		Timer timer_H;			 				/**< total time of Hessian computation */
		Timer timer_fs; 						/**< total time of manipulation with function values during iterations */
		Timer timer_integrate;	 				/**< total time of integration */

		EntropyData<VectorBase> *entropydata; 	/**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments from data, size K*Km */
		GeneralVector<VectorBase> *integrals; /**< vector of computed integrals, size K*(2*Km+1) */
		GeneralVector<VectorBase> *x_power; /**< global temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< global temp vector for storing power of x * gamma_k */

		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();

		int debugmode;				/**< basic debug mode schema [0/1/2] */
		bool debug_print_it;		/**< print simple info about outer iterations */
		bool debug_print_vectors;	/**< print content of vectors during iterations */
		bool debug_print_scalars;	/**< print values of computed scalars during iterations */ 
		bool debug_print_Axb;		/**< print info about inner Ax=b solver every outer Newton iteration */

		bool monitor;				/**< export the descend into .m file */

		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		void prepare_entropyintegration(int integration_type, double integration_eps);

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		/* Ax=b stuff (for problem of size Km) */
		GeneralVector<VectorBase> *g; 		/**< local gradient, size Km */
		GeneralVector<VectorBase> *delta;	/**< vetor used in Newton method, size Km */

	public:

		EntropySolverNewton();
		EntropySolverNewton(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolverNewton();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printcontent(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		EntropyData<VectorBase> *get_data() const;

		void compute_moments();
		void compute_residuum(GeneralVector<VectorBase> *residuum) const;
		
		ExternalContent *get_externalcontent() const;

		int get_xdim() const;
		int get_K() const;
		int get_Km() const;
		int get_number_of_moments() const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

template<class VectorBase>
void EntropySolverNewton<VectorBase>::prepare_entropyintegration(int integration_type, double integration_eps) {
	LOG_FUNC_BEGIN

	/* auto */
	if(integration_type == 0){
		integration_type = 1;
	}

	/* dlib */
	if(integration_type == 1){
		/* can be used only for xdim=1 */
		//TODO: if(this->xdim != 1) throw error 
		this->entropyintegration = new EntropyIntegrationDlib<VectorBase>(this->entropydata, integration_eps);
	}

	/* cuba */
	if(integration_type == 2){
		this->entropyintegration = new EntropyIntegrationCuba<VectorBase>(this->entropydata, integration_eps);
	}

	LOG_FUNC_END
}


template<class VectorBase>
void EntropySolverNewton<VectorBase>::set_settings_from_console() {
	LOG_FUNC_BEGIN

	consoleArg.set_option_value("entropysolvernewton_maxit", &this->maxit, ENTROPYSOLVERNEWTON_DEFAULT_MAXIT);
	consoleArg.set_option_value("entropysolvernewton_maxit_Axb", &this->maxit_Axb, ENTROPYSOLVERNEWTON_DEFAULT_MAXIT_AXB);
	consoleArg.set_option_value("entropysolvernewton_eps", &this->eps, ENTROPYSOLVERNEWTON_DEFAULT_EPS);
	consoleArg.set_option_value("entropysolvernewton_eps_Axb", &this->eps_Axb, ENTROPYSOLVERNEWTON_DEFAULT_EPS_AXB);
	consoleArg.set_option_value("entropysolvernewton_newton_coeff", &this->newton_coeff, ENTROPYSOLVERNEWTON_DEFAULT_NEWTON_COEFF);
	
	int integration_type;
	double integration_eps;
	consoleArg.set_option_value("entropysolvernewton_integration_eps", &integration_eps, ENTROPYSOLVERNEWTON_DEFAULT_INTEGRATION_EPS);
	consoleArg.set_option_value("entropysolvernewton_integration_type", &integration_type, ENTROPYSOLVERNEWTON_DEFAULT_INTEGRATION_TYPE);
	prepare_entropyintegration(integration_type, integration_eps);
	
	/* based on integration type allocate EntropyIntegration */
	consoleArg.set_option_value("entropysolvernewton_monitor", &this->monitor, ENTROPYSOLVERNEWTON_MONITOR);	

	/* set debug mode */
	consoleArg.set_option_value("entropysolvernewton_debugmode", &this->debugmode, ENTROPYSOLVERNEWTON_DEFAULT_DEBUGMODE);

	
	debug_print_vectors = false;
	debug_print_scalars = false; 
	debug_print_it = false; 

	if(debugmode == 1){
		debug_print_it = true;
	}

	if(debugmode == 2){
		debug_print_it = true;
		debug_print_scalars = true;
		debug_print_Axb = true;
	}

	consoleArg.set_option_value("entropysolvernewton_debug_print_it",		&debug_print_it, 		debug_print_it);
	consoleArg.set_option_value("entropysolvernewton_debug_print_vectors", 	&debug_print_vectors,	false);
	consoleArg.set_option_value("entropysolvernewton_debug_print_scalars", 	&debug_print_scalars, 	debug_print_scalars);
	consoleArg.set_option_value("entropysolvernewton_debug_print_Axb", 		&debug_print_Axb, 		debug_print_Axb);

	LOG_FUNC_END
}

/* prepare temp_vectors */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* constructor */
template<class VectorBase>
EntropySolverNewton<VectorBase>::EntropySolverNewton(){
	LOG_FUNC_BEGIN
	
	entropydata = NULL;

	/* initial values */
	this->it_sums = new int [1];
	this->it_lasts = new int [1];
	this->itAxb_sums = new int [1];
	this->itAxb_lasts = new int [1];
	this->fxs = new double [1];
	this->gnorms = new double [1];

	set_value_array(1, this->it_sums, 0);
	set_value_array(1, this->it_lasts, 0);
	set_value_array(1, this->itAxb_sums, 0);
	set_value_array(1, this->itAxb_lasts, 0);
	set_value_array(1, this->fxs, std::numeric_limits<double>::max());
	set_value_array(1, this->gnorms, std::numeric_limits<double>::max());

	this->number_of_moments = 0;
	this->xdim = 1;
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_compute_moments.restart();	
	this->timer_solve.restart();	
	this->timer_Axb.restart();	
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_H.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	LOG_FUNC_END
}

template<class VectorBase>
EntropySolverNewton<VectorBase>::EntropySolverNewton(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;
	int K = entropydata->get_K();

	/* initial values */
	this->it_sums = new int [K];
	this->it_lasts = new int [K];
	this->itAxb_sums = new int [K];
	this->itAxb_lasts = new int [K];
	this->fxs = new double [K];
	this->gnorms = new double [K];

	set_value_array(K, this->it_sums, 0);
	set_value_array(K, this->it_lasts, 0);
	set_value_array(K, this->itAxb_sums, 0);
	set_value_array(K, this->itAxb_lasts, 0);
	set_value_array(K, this->fxs, std::numeric_limits<double>::max());
	set_value_array(K, this->gnorms, std::numeric_limits<double>::max());

	this->number_of_moments = EntropyData<VectorBase>::compute_number_of_moments(get_xdim(), get_Km());
	this->xdim = entropydata->get_xdim();

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_compute_moments.restart();	
	this->timer_solve.restart();
	this->timer_Axb.restart();
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_H.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	allocate_temp_vectors();
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolverNewton<VectorBase>::~EntropySolverNewton(){
	LOG_FUNC_BEGIN

	free(this->it_sums);
	free(this->it_lasts);
	free(this->itAxb_sums);
	free(this->itAxb_lasts);
	free(this->fxs);
	free(this->gnorms);

	/* free temp vectors */
	free_temp_vectors();

	/* free tool for integration */
	if(this->entropyintegration){
		free(this->entropyintegration);
	}

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit             : " << this->maxit << std::endl;
	output <<  " - maxit_Axb         : " << this->maxit_Axb << std::endl;
	output <<  " - eps               : " << this->eps << std::endl;
	output <<  " - eps_Axb           : " << this->eps_Axb << std::endl;
	output <<  " - newton_coeff      : " << this->newton_coeff << std::endl;
	output <<  " - integrationtype   : " << this->entropyintegration->get_name() << std::endl;
	output.push();
	this->entropyintegration->print(output);
	output.pop();
	
	output <<  " - xdim              : " << this->get_xdim() << std::endl;
	output <<  " - K                 : " << this->get_K() << std::endl;
	output <<  " - Km                : " << this->get_Km() << std::endl;
	output <<  " - number_of_moments : " << this->get_number_of_moments() << std::endl;
	
	output <<  " - debugmode         : " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		coutMaster.push();
		entropydata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit             : " << this->maxit << std::endl;
	output_global <<  " - maxit_Axb         : " << this->maxit_Axb << std::endl;
	output_global <<  " - eps               : " << this->eps << std::endl;
	output_global <<  " - eps_Axb           : " << this->eps_Axb << std::endl;
	output_global <<  " - newton_coeff      : " << this->newton_coeff << std::endl;
	output_global <<  " - integrationtype   : " << this->entropyintegration->get_name() << std::endl;
	output_global.push();
	this->entropyintegration->print(output_global, output_local);
	output_global.pop();

	output_global <<  " - xdim              : " << this->get_xdim() << std::endl;
	output_global <<  " - K                 : " << this->get_K() << std::endl;
	output_global <<  " - Km                : " << this->get_Km() << std::endl;
	output_global <<  " - number_of_moments : " << this->get_number_of_moments() << std::endl;

	output_global <<  " - debugmode         : " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		output_global << "- data:" << std::endl;
		coutMaster.push();
		entropydata->print(output_global,output_local);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	int K = entropydata->get_K();
	
	output <<  this->get_name() << std::endl;
	output <<  " - it: " << std::setw(6) << print_array(this->it_lasts, K) << ", ";
	output <<  "fx: " << std::setw(10) << print_array(this->fxs, K) << ", ";	
	output <<  "norm(g): " << std::setw(10) << print_array(this->gnorms, K) << ", ";
	output <<  "used memory: " << std::setw(6) << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	int K = entropydata->get_K();
	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - fx          : " << std::setw(25) << print_array(this->fxs, K) << std::endl;
	output <<  "      - norm(g)     : " << std::setw(25) << print_array(this->gnorms, K) << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(entropydata){
		output << "- data:" << std::endl;
		coutMaster.push();
		entropydata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	int K = entropydata->get_K();

	output <<  this->get_name() << std::endl;
	output <<  " - it all        = " << print_array(this->it_sums,K) << std::endl;
	output <<  " - itAxb all     = " << print_array(this->itAxb_sums,K) << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_moments    = " << this->timer_compute_moments.get_value_sum() << std::endl;
	output <<  "  - t_solve      = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_Axb        = " << this->timer_Axb.get_value_sum() << std::endl;
	output <<  "  - t_update     = " << this->timer_update.get_value_sum() << std::endl;
	output <<  "  - t_g          = " << this->timer_g.get_value_sum() << std::endl;
	output <<  "  - t_H          = " << this->timer_H.get_value_sum() << std::endl;
	output <<  "  - t_fs         = " << this->timer_fs.get_value_sum() << std::endl;
	output <<  "  - t_integrate  = " << this->timer_integrate.get_value_sum() << std::endl;
	output <<  "  - t_other      = " << this->timer_solve.get_value_sum() - (this->timer_integrate.get_value_sum() + this->timer_H.get_value_sum() + this->timer_update.get_value_sum() + this->timer_g.get_value_sum() + this->timer_fs.get_value_sum() + this->timer_Axb.get_value_sum()) << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolverNewton<VectorBase>::get_name() const {
	std::string return_value = "EntropySolverNewton<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
EntropyData<VectorBase> *EntropySolverNewton<VectorBase>::get_data() const {
	return entropydata;
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::compute_moments() {
	LOG_FUNC_BEGIN

	entropydata->compute_moments(this->moments);

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::compute_residuum(GeneralVector<VectorBase> *residuum) const {
	LOG_FUNC_BEGIN

	entropydata->compute_residuum(residuum, this->integrals);
		
	LOG_FUNC_END
}

template<class VectorBase>
int EntropySolverNewton<VectorBase>::get_xdim() const {
	return entropydata->get_xdim();
}

template<class VectorBase>
int EntropySolverNewton<VectorBase>::get_K() const {
	return entropydata->get_K();
}

template<class VectorBase>
int EntropySolverNewton<VectorBase>::get_Km() const {
	return entropydata->get_Km();
}

template<class VectorBase>
int EntropySolverNewton<VectorBase>::get_number_of_moments() const {
	return entropydata->get_number_of_moments();
}


}
} /* end namespace */

#endif
