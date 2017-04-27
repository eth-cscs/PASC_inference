/** @file entropysolverspg.h
 *  @brief Solver which solves problem with integrals using SPG algorithm
 *  
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERSPG_H
#define	PASC_ENTROPYSOLVERSPG_H

#include "general/solver/generalsolver.h"
#include "general/data/entropydata.h"
#include "general/solver/spg_fs.h"

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

#define ENTROPYSOLVERSPG_DEFAULT_MAXIT 1000
#define ENTROPYSOLVERSPG_DEFAULT_MAXIT_GLL 100
#define ENTROPYSOLVERSPG_DEFAULT_EPS 1e-6
#define ENTROPYSOLVERSPG_DEFAULT_DEBUGMODE 0

#define ENTROPYSOLVERSPG_DEFAULT_M 20
#define ENTROPYSOLVERSPG_DEFAULT_GAMMA 0.9
#define ENTROPYSOLVERSPG_DEFAULT_SIGMA1 0.0
#define ENTROPYSOLVERSPG_DEFAULT_SIGMA2 1.0
#define ENTROPYSOLVERSPG_DEFAULT_ALPHAINIT 0.1

#define ENTROPYSOLVERSPG_MONITOR false

namespace pascinference {
namespace solver {

/** \class EntropySolverSPG
 *  \brief Solver which solves problem with integrals using SPG algorithm
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolverSPG: public GeneralSolver {
	protected:
		double *fxs;					/**< function values in clusters */
		double *gnorms;					/**< norm of gradient in clusters (stopping criteria) */
	
		int maxit_gll;
		int *it_sums;					/**< sums of all iterations for each cluster */
		int *it_lasts;					/**< number of iterations from last solve() call for each cluster */
		int *itgll_sums;				/**< sums of all gll iterations for each cluster */
		int *itgll_lasts;				/**< sums of all gll iterations in this outer iteration */
	
		Timer timer_compute_moments;	/**< time for computing moments from data */
		Timer timer_solve; 				/**< total solution time of SPG algorithm */
		Timer timer_update; 			/**< total time of vector updates */
		Timer timer_g;			 		/**< total time of gradient computation */
		Timer timer_dot;		 		/**< total time of dot computation */
		Timer timer_fs; 				/**< total time of manipulation with fs vector during iterations */
		Timer timer_integrate;	 		/**< total time of integration */

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

		bool monitor;				/**< export the descend into .m file */

		int m;						/**< size of SPG_fs */
		double gamma;				/**< parameter of Armijo condition */
		double sigma1;				/**< to enforce progress */
		double sigma2;				/**< to enforce progress */
		double alphainit;			/** initial step-size */

		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		/* SPG stuff (for problem of size Km) */
		GeneralVector<VectorBase> *g; 		/**< local gradient, size Km */
		GeneralVector<VectorBase> *y; 		/**< g_old, g-g_old, size Km */
		GeneralVector<VectorBase> *s; 		/**< x_old, x-x_old, size Km */

//TODO: !!!!
#ifdef USE_PETSC
		void compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);					/**< g = nabla_lambda f(lambda, moments) */
		double compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);		/**< compute function value from already computed integrals and moments */
		void compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all);				/**< compute integrals int x^{0,..,Km} exp(-dot(lambda,x^{1,..,Km})) */
#endif

	public:

		EntropySolverSPG();
		EntropySolverSPG(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolverSPG();

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
namespace pascinference {
namespace solver {

template<class VectorBase>
void EntropySolverSPG<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("entropysolverspg_maxit", &this->maxit, ENTROPYSOLVERSPG_DEFAULT_MAXIT);
	consoleArg.set_option_value("entropysolverspg_maxit_gll", &this->maxit_gll, ENTROPYSOLVERSPG_DEFAULT_MAXIT_GLL);
	consoleArg.set_option_value("entropysolverspg_eps", &this->eps, ENTROPYSOLVERSPG_DEFAULT_EPS);
	
	consoleArg.set_option_value("entropysolverspg_m", &this->m, ENTROPYSOLVERSPG_DEFAULT_M);	
	consoleArg.set_option_value("entropysolverspg_gamma", &this->gamma, ENTROPYSOLVERSPG_DEFAULT_GAMMA);	
	consoleArg.set_option_value("entropysolverspg_sigma1", &this->sigma1, ENTROPYSOLVERSPG_DEFAULT_SIGMA1);	
	consoleArg.set_option_value("entropysolverspg_sigma2", &this->sigma2, ENTROPYSOLVERSPG_DEFAULT_SIGMA2);	
	consoleArg.set_option_value("entropysolverspg_alphainit", &this->alphainit, ENTROPYSOLVERSPG_DEFAULT_ALPHAINIT);	

	consoleArg.set_option_value("entropysolverspg_monitor", &this->monitor, ENTROPYSOLVERSPG_MONITOR);	

	/* set debug mode */
	consoleArg.set_option_value("entropysolverspg_debugmode", &this->debugmode, ENTROPYSOLVERSPG_DEFAULT_DEBUGMODE);

	
	debug_print_vectors = false;
	debug_print_scalars = false; 
	debug_print_it = false; 

	if(debugmode == 1){
		debug_print_it = true;
	}

	if(debugmode == 2){
		debug_print_it = true;
		debug_print_scalars = true;
	}

	consoleArg.set_option_value("entropysolverspg_debug_print_it",		&debug_print_it, 		debug_print_it);
	consoleArg.set_option_value("entropysolverspg_debug_print_vectors", &debug_print_vectors,	false);
	consoleArg.set_option_value("entropysolverspg_debug_print_scalars", &debug_print_scalars, 	debug_print_scalars);

}

/* prepare temp_vectors */
template<class VectorBase>
void EntropySolverSPG<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void EntropySolverSPG<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


/* constructor */
template<class VectorBase>
EntropySolverSPG<VectorBase>::EntropySolverSPG(){
	LOG_FUNC_BEGIN
	
	entropydata = NULL;

	/* initial values */
	this->it_sums = new int [1];
	this->it_lasts = new int [1];
	this->itgll_sums = new int [1];
	this->itgll_lasts = new int [1];
	this->fxs = new double [1];
	this->gnorms = new double [1];

	set_value_array(1, this->it_sums, 0);
	set_value_array(1, this->it_lasts, 0);
	set_value_array(1, this->itgll_sums, 0);
	set_value_array(1, this->itgll_lasts, 0);
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
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_dot.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	LOG_FUNC_END
}

template<class VectorBase>
EntropySolverSPG<VectorBase>::EntropySolverSPG(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;
	int K = entropydata->get_K();

	/* initial values */
	this->it_sums = new int [K];
	this->it_lasts = new int [K];
	this->itgll_sums = new int [K];
	this->itgll_lasts = new int [K];
	this->fxs = new double [K];
	this->gnorms = new double [K];

	set_value_array(K, this->it_sums, 0);
	set_value_array(K, this->it_lasts, 0);
	set_value_array(K, this->itgll_sums, 0);
	set_value_array(K, this->itgll_lasts, 0);
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
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_dot.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	allocate_temp_vectors();
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolverSPG<VectorBase>::~EntropySolverSPG(){
	LOG_FUNC_BEGIN

	free(this->it_sums);
	free(this->it_lasts);
	free(this->itgll_sums);
	free(this->itgll_lasts);	
	free(this->fxs);
	free(this->gnorms);

	/* free temp vectors */
	free_temp_vectors();

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolverSPG<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - maxit_gll:  " << this->maxit_gll << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	output <<  " - m:          " << m << std::endl;
	output <<  " - gamma:      " << gamma << std::endl;
	output <<  " - sigma1:     " << sigma1 << std::endl;
	output <<  " - sigma2:     " << sigma2 << std::endl;
	output <<  " - alphainit:  " << alphainit << std::endl;

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
void EntropySolverSPG<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << std::endl;
	output_global <<  " - maxit_gll:  " << this->maxit_gll << std::endl;
	output_global <<  " - eps:        " << this->eps << std::endl;
	output_global <<  " - debugmode: " << this->debugmode << std::endl;

	output_local <<  " - m:          " << m << std::endl;
	output_local <<  " - gamma:      " << gamma << std::endl;
	output_local <<  " - sigma1:     " << sigma1 << std::endl;
	output_local <<  " - sigma2:     " << sigma2 << std::endl;
	output_local <<  " - alphainit:  " << alphainit << std::endl;

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
void EntropySolverSPG<VectorBase>::printstatus(ConsoleOutput &output) const {
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
void EntropySolverSPG<VectorBase>::printstatus(std::ostringstream &output) const {
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
void EntropySolverSPG<VectorBase>::printcontent(ConsoleOutput &output) const {
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
void EntropySolverSPG<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolverSPG<VectorBase>::get_name() const {
	std::string return_value = "EntropySolverSPG<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
void EntropySolverSPG<VectorBase>::solve() {
	LOG_FUNC_BEGIN
	
	//TODO

	LOG_FUNC_END
}	

template<class VectorBase>
void EntropySolverSPG<VectorBase>::compute_moments_data() {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

#ifdef USE_DLIB
template<class VectorBase>
double EntropySolverSPG<VectorBase>::gg(double y, int order, column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    
    return pow(y,order)*(exp(-trans(LM)*z));
}
#endif

template<class VectorBase>
void EntropySolverSPG<VectorBase>::compute_residuum(GeneralVector<VectorBase> *residuum) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


}
} /* end namespace */

#endif
