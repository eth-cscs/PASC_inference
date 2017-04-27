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
#include "general/algebra/integration/entropyintegration.h"

/* include Dlib stuff */
#ifdef USE_DLIB
	#include "dlib/matrix.h"
	#include "dlib/numeric_constants.h"
	#include "dlib/numerical_integration.h"
	#include "dlib/optimization.h"

	/* Dlib column vector */
	typedef dlib::matrix<double,0,1> column_vector;
#endif


#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT 1000
#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT_CG 100
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS 1e-6
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS_CG 1e-6
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
	protected:
		/** @brief type of numerical integration
		 * 
		 */
		typedef enum { 
			INTEGRATION_AUTO=0,					/**< choose automatic solver */
			INTEGRATION_DLIB=1,					/**< use Dlib library to compute integrals */
			INTEGRATION_MC=2					/**< use Monte Carlo integration method */
		} IntegrationType;

		/** @brief return name of integration solver in string format
		 */
		std::string print_integrationtype(IntegrationType integrationtype_in) const;
	
		double *fxs;							/**< function values in clusters */
		double *gnorms;							/**< norm of gradient in clusters (stopping criteria) */
	
		int maxit_ksp;
		double eps_ksp;
		int *it_sums;							/**< sums of all iterations for each cluster */
		int *it_lasts;							/**< number of iterations from last solve() call for each cluster */
		int *itksp_sums;						/**< sums of all cg iterations for each cluster */
		int *itksp_lasts;						/**< sums of all cg iterations in this outer iteration */
		double newton_coeff;					/**< newton step-size coefficient x_{k+1} = x_k + coeff*delta */
		IntegrationType integrationtype;	 	/**< the type of numerical integration */
		EntropyIntegration<VectorBase> *entropyintegration;	/**< instance of integration tool */
	
		Timer timer_compute_moments;			/**< time for computing moments from data */
		Timer timer_solve; 						/**< total solution time of Newton algorithm */
		Timer timer_ksp; 						/**< total solution time of KSP algorithm */
		Timer timer_update; 					/**< total time of vector updates */
		Timer timer_g;			 				/**< total time of gradient computation */
		Timer timer_H;			 				/**< total time of Hessian computation */
		Timer timer_fs; 						/**< total time of manipulation with function values during iterations */
		Timer timer_integrate;	 				/**< total time of integration */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments_data; /**< vector of computed moments from data, size K*Km */
		GeneralVector<VectorBase> *integrals; /**< vector of computed integrals, size K*(Km+1) */
		GeneralVector<VectorBase> *x_power; /**< global temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< global temp vector for storing power of x * gamma_k */

		/* functions for Dlib */
		#ifdef USE_DLIB
			static double gg(double y, int order, column_vector& LM);
		#endif

		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();
		
		int debugmode;				/**< basic debug mode schema [0/1/2] */
		bool debug_print_it;		/**< print simple info about outer iterations */
		bool debug_print_vectors;	/**< print content of vectors during iterations */
		bool debug_print_scalars;	/**< print values of computed scalars during iterations */ 
		bool debug_print_ksp;		/**< print info about inner KSP every outer Newton iteration */

		bool monitor;				/**< export the descend into .m file */

		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		/* KSP stuff (for problem of size Km) */
		GeneralVector<VectorBase> *g; 		/**< local gradient, size Km */
		GeneralVector<VectorBase> *delta;	/**< vetor used in Newton method, size Km */

// TODO: !!!!
#ifdef USE_PETSC
		Mat H_petsc;						/**< Hessian matrix */
		KSP ksp;							/**< linear solver context */
		PC pc;           					/**< preconditioner context **/

		void compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);					/**< g = nabla_lambda f(lambda, moments) */
		void compute_hessian(Vec &integrals_Vec);													/**< Hessian matrix */
		double compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);		/**< compute function value from already computed integrals and moments */
		void compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all);				/**< compute integrals int x^{0,..,Km} exp(-dot(lambda,x^{1,..,Km})) */

#endif		

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

		void compute_moments_data();
		
		void compute_residuum(GeneralVector<VectorBase> *residuum) const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

template<class VectorBase>
std::string EntropySolverNewton<VectorBase>::print_integrationtype(IntegrationType integrationtype_in) const{
	std::string return_value = "?";
	switch(integrationtype_in){
		case(INTEGRATION_AUTO): return_value = "AUTO"; break;
		case(INTEGRATION_DLIB): return_value = "DLIB"; break;
		case(INTEGRATION_MC):   return_value = "Monte Carlo"; break;
	}
	return return_value;
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("entropysolvernewton_maxit", &this->maxit, ENTROPYSOLVERNEWTON_DEFAULT_MAXIT);
	consoleArg.set_option_value("entropysolvernewton_maxit_ksp", &this->maxit_ksp, ENTROPYSOLVERNEWTON_DEFAULT_MAXIT_CG);
	consoleArg.set_option_value("entropysolvernewton_eps", &this->eps, ENTROPYSOLVERNEWTON_DEFAULT_EPS);
	consoleArg.set_option_value("entropysolvernewton_eps_ksp", &this->eps_ksp, ENTROPYSOLVERNEWTON_DEFAULT_EPS_CG);
	consoleArg.set_option_value("entropysolvernewton_newton_coeff", &this->newton_coeff, ENTROPYSOLVERNEWTON_DEFAULT_NEWTON_COEFF);
	
	int integrationtype_int;
	consoleArg.set_option_value("entropysolvernewton_integrationtype", &integrationtype_int, INTEGRATION_AUTO);
	this->integrationtype = static_cast<IntegrationType>(integrationtype_int);
	
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
		debug_print_ksp = true;
	}

	consoleArg.set_option_value("entropysolvernewton_debug_print_it",		&debug_print_it, 		debug_print_it);
	consoleArg.set_option_value("entropysolvernewton_debug_print_vectors", 	&debug_print_vectors,	false);
	consoleArg.set_option_value("entropysolvernewton_debug_print_scalars", 	&debug_print_scalars, 	debug_print_scalars);
	consoleArg.set_option_value("entropysolvernewton_debug_print_ksp", 		&debug_print_ksp, 		debug_print_ksp);

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
	this->itksp_sums = new int [1];
	this->itksp_lasts = new int [1];
	this->fxs = new double [1];
	this->gnorms = new double [1];

	set_value_array(1, this->it_sums, 0);
	set_value_array(1, this->it_lasts, 0);
	set_value_array(1, this->itksp_sums, 0);
	set_value_array(1, this->itksp_lasts, 0);
	set_value_array(1, this->fxs, std::numeric_limits<double>::max());
	set_value_array(1, this->gnorms, std::numeric_limits<double>::max());
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_compute_moments.restart();	
	this->timer_solve.restart();	
	this->timer_ksp.restart();	
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_H.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	this->entropyintegration = NULL;

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
	this->itksp_sums = new int [K];
	this->itksp_lasts = new int [K];
	this->fxs = new double [K];
	this->gnorms = new double [K];

	set_value_array(K, this->it_sums, 0);
	set_value_array(K, this->it_lasts, 0);
	set_value_array(K, this->itksp_sums, 0);
	set_value_array(K, this->itksp_lasts, 0);
	set_value_array(K, this->fxs, std::numeric_limits<double>::max());
	set_value_array(K, this->gnorms, std::numeric_limits<double>::max());

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_compute_moments.restart();	
	this->timer_solve.restart();
	this->timer_ksp.restart();
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_H.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	allocate_temp_vectors();

	this->entropyintegration = NULL; /* not now, the integrator will be initialized with first integration call */
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolverNewton<VectorBase>::~EntropySolverNewton(){
	LOG_FUNC_BEGIN

	free(this->it_sums);
	free(this->it_lasts);
	free(this->itksp_sums);
	free(this->itksp_lasts);
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
	output <<  " - maxit            : " << this->maxit << std::endl;
	output <<  " - maxit_ksp        : " << this->maxit_ksp << std::endl;
	output <<  " - eps              : " << this->eps << std::endl;
	output <<  " - eps_ksp          : " << this->eps_ksp << std::endl;
	output <<  " - newton_coeff     : " << this->newton_coeff << std::endl;
	output <<  " - integrationtype  : " << print_integrationtype(this->integrationtype) << std::endl;
	
	output <<  " - debugmode        : " << this->debugmode << std::endl;

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
	output_global <<  " - maxit            : " << this->maxit << std::endl;
	output_global <<  " - maxit_ksp        : " << this->maxit_ksp << std::endl;
	output_global <<  " - eps              : " << this->eps << std::endl;
	output_global <<  " - eps_ksp          : " << this->eps_ksp << std::endl;
	output_global <<  " - newton_coeff     : " << this->newton_coeff << std::endl;
	output_global <<  " - integrationtype  : " << print_integrationtype(this->integrationtype) << std::endl;

	output_global <<  " - debugmode    : " << this->debugmode << std::endl;

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
	output <<  " - itksp all     = " << print_array(this->itksp_sums,K) << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_moments    = " << this->timer_compute_moments.get_value_sum() << std::endl;
	output <<  "  - t_solve      = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_ksp        = " << this->timer_ksp.get_value_sum() << std::endl;
	output <<  "  - t_update     = " << this->timer_update.get_value_sum() << std::endl;
	output <<  "  - t_g          = " << this->timer_g.get_value_sum() << std::endl;
	output <<  "  - t_H          = " << this->timer_H.get_value_sum() << std::endl;
	output <<  "  - t_fs         = " << this->timer_fs.get_value_sum() << std::endl;
	output <<  "  - t_integrate  = " << this->timer_integrate.get_value_sum() << std::endl;
	output <<  "  - t_other      = " << this->timer_solve.get_value_sum() - (this->timer_integrate.get_value_sum() + this->timer_H.get_value_sum() + this->timer_update.get_value_sum() + this->timer_g.get_value_sum() + this->timer_fs.get_value_sum() + this->timer_ksp.get_value_sum()) << std::endl;

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
void EntropySolverNewton<VectorBase>::compute_moments_data() {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverNewton<VectorBase>::compute_residuum(GeneralVector<VectorBase> *residuum) const {
	LOG_FUNC_BEGIN

	//TODO
		
	LOG_FUNC_END
}

#ifdef USE_DLIB
template<class VectorBase>
double EntropySolverNewton<VectorBase>::gg(double y, int order, column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    
    return pow(y,order)*(exp(-trans(LM)*z));
}
#endif


}
} /* end namespace */

#endif
