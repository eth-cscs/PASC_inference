/** @file entropysolvernewton.h
 *  @brief Solver which solves problem with integrals using Newton method
 *  
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERNEWTON_H
#define	PASC_ENTROPYSOLVERNEWTON_H

#include "pascinference.h"
#include "data/entropydata.h"

/* include integration algorithms */
#include "algebra/integration/entropyintegration.h"

/* petsc stuff */
#include <petscksp.h>

/* include Dlib stuff */
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT 1000
#define ENTROPYSOLVERNEWTON_DEFAULT_MAXIT_CG 100
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS 1e-6
#define ENTROPYSOLVERNEWTON_DEFAULT_EPS_CG 1e-6
#define ENTROPYSOLVERNEWTON_DEFAULT_NEWTON_COEFF 0.9
#define ENTROPYSOLVERNEWTON_DEFAULT_DEBUGMODE 0

#define ENTROPYSOLVERNEWTON_MONITOR false

/* Dlib column vector */
typedef dlib::matrix<double,0,1> column_vector;

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
		EntropyIntegration *entropyintegration;	/**< instance of integration tool */
	
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
		GeneralVector<VectorBase> *x_gammak; /**< global temp vector for storing x*gamma_k */
		GeneralVector<VectorBase> *x_gammak_power; /**< global temp vector for storing power of x * gamma_k */

		/* functions for Dlib */
		static double gg(double y, int order, column_vector& LM);

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
		Mat H_petsc;						/**< Hessian matrix */
		KSP ksp;							/**< linear solver context */
		PC pc;           					/**< preconditioner context **/
		

		void compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);					/**< g = nabla_lambda f(lambda, moments) */
		void compute_hessian(Vec &integrals_Vec);													/**< Hessian matrix */
		double compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);		/**< compute function value from already computed integrals and moments */
		void compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all);				/**< compute integrals int x^{0,..,Km} exp(-dot(lambda,x^{1,..,Km})) */

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

	/* prepare auxiliary vectors */
	x_gammak = new GeneralVector<PetscVector>(*entropydata->get_x());
	x_gammak_power = new GeneralVector<PetscVector>(*entropydata->get_x());
	
	/* create aux vector for the computation of moments and integrals */
	Vec moments_Vec;
	Vec integrals_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&moments_Vec) );
	TRYCXX( VecCreate(PETSC_COMM_SELF,&integrals_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(moments_Vec, VECSEQCUDA));
		TRYCXX(VecSetType(integrals_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(moments_Vec, VECSEQ));
		TRYCXX(VecSetType(integrals_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*entropydata->get_Km(),PETSC_DECIDE) );
	TRYCXX( VecSetSizes(integrals_Vec,entropydata->get_K()*(2*entropydata->get_Km()+1),PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	TRYCXX( VecSetFromOptions(integrals_Vec) );
	this->moments_data = new GeneralVector<PetscVector>(moments_Vec);
	this->integrals = new GeneralVector<PetscVector>(integrals_Vec);

	/* SPG stuff - create vector to solve k-th problem (of size Km) */
	Vec g_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&g_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(g_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(g_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(g_Vec,entropydata->get_Km(),entropydata->get_Km()) );
	TRYCXX( VecSetFromOptions(g_Vec) );
	this->g = new GeneralVector<PetscVector>(g_Vec);
	this->delta = new GeneralVector<PetscVector>(*g);


	/* create Hessian matrix */
	TRYCXX( MatCreate(PETSC_COMM_SELF, &H_petsc) );
	TRYCXX( MatSetSizes(H_petsc,entropydata->get_Km(),entropydata->get_Km(),PETSC_DECIDE,PETSC_DECIDE) );
	#ifdef USE_CUDA
		TRYCXX( MatSetType(H_petsc,MATSEQAIJCUSPARSE) ); 
	#else
		TRYCXX( MatSetType(H_petsc,MATSEQAIJ) ); 
	#endif
	TRYCXX( MatSetFromOptions(H_petsc) );
	TRYCXX( MatMPIAIJSetPreallocation(H_petsc,entropydata->get_Km(),NULL,entropydata->get_Km()-1,NULL) ); 
	TRYCXX( MatSeqAIJSetPreallocation(H_petsc,entropydata->get_Km(),NULL) );

	TRYCXX( MatAssemblyBegin(H_petsc,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(H_petsc,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)H_petsc,"Hessian matrix") );

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void EntropySolverNewton<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(x_gammak);
	free(x_gammak_power);
	free(moments_data);
	free(integrals);

	free(g);
	free(delta);

	TRYCXX( MatDestroy(&H_petsc) );

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

	this->timer_compute_moments.start();
	 this->compute_moments_data();
	this->timer_compute_moments.stop();
	
//	coutMaster << "Moments: " << *moments << std::endl;

	this->timer_solve.start(); 
	int it; /* actual number of iterations */
	int itksp, itksp_one_newton_iteration; /* number of all KSP iterations, number of KSP iterations in one newton iteration */
	
	/* get dimensions */
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* get PETSc vecs - LOCAL */
	Vec moments_Vec = moments_data->get_vector();
	Vec x_Vec = entropydata->get_lambda()->get_vector(); /* x:=lambda is unknowns */
	Vec g_Vec = g->get_vector();
	Vec delta_Vec = delta->get_vector();
	Vec integrals_Vec = integrals->get_vector();

	/* stuff for clusters (subvectors) - LOCAL */
	IS k_is;
	IS integralsk_is;
	Vec xk_Vec;
	Vec momentsk_Vec;
	Vec integralsk_Vec;

	double gnorm, fx; /* norm(g),f(x) */
	double deltanorm; /* norm(delta) - for debug */

	/* postprocess */
	Vec g_inner_Vec;
	double gnorm_inner;

	/* through all clusters */
	for(int k = 0; k < K; k++){
		
		coutMaster << "k=" << k << std::endl;
		
		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km, k*Km, 1, &k_is) ); /* Theta is LOCAL ! */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km+1, k*(Km+1), 1, &integralsk_is) ); 
	
		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(x_Vec, k_is, &xk_Vec) );
		TRYCXX( VecGetSubVector(moments_Vec, k_is, &momentsk_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* -------------- Newton algorithm (for k-th problem) ------------ */
		it = 0;
		itksp = 0;
		
		/* compute integrals and gradient */
		this->timer_integrate.start();
		 compute_integrals(integralsk_Vec, xk_Vec, true);
		this->timer_integrate.stop();
		this->timer_g.start();
		 compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_g.stop();
		this->timer_fs.start();
		 fx = compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_fs.stop();

		while(it < this->maxit){
			/* compute stopping criteria - norm of gradient */
			TRYCXX( VecNorm(g_Vec, NORM_2, &gnorm) );
			if(gnorm < this->eps){
				break;
			}

			/* prepare Hessian matrix */
			this->timer_H.start();
			 compute_hessian(integralsk_Vec);
			this->timer_H.stop();
			
			/* ------ KSP Solver ----- */
			/* for solving H*delta=-g */
			TRYCXX( VecScale(g_Vec,-1.0) );
			
			/* Create linear solver context */
			TRYCXX( KSPCreate(PETSC_COMM_SELF,&ksp) ); // TODO: really every iteration?
			
			/* Set operators. Here the matrix that defines the linear system */
			/* also serves as the preconditioning matrix. */
			TRYCXX( KSPSetOperators(ksp,H_petsc,H_petsc) ); 
			
			/* precondition - maybe we can try LU factorization - then system will be solved in one iteration */
			TRYCXX( KSPGetPC(ksp,&pc) );
			TRYCXX( PCSetType(pc,PCJACOBI) );
			
			/* set stopping criteria */
			TRYCXX( KSPSetTolerances(ksp,this->eps_ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT) ); // TODO: look onto all arguments

			/* some funny stuff can be loaded from petsc console parameters */
			TRYCXX( KSPSetFromOptions(ksp) );
			
			/* I think that these thing will use actual values in delta as initial approximation */
			TRYCXX( KSPSetInitialGuessNonzero(ksp,PETSC_TRUE) );
			
			/* Solve linear system */
			this->timer_ksp.start();
			 TRYCXX( KSPSolve(ksp,g_Vec,delta_Vec) );
			this->timer_ksp.stop();
			
			/* get some informations about solution process */
			TRYCXX( KSPGetIterationNumber(ksp,&itksp_one_newton_iteration) );
			itksp += itksp_one_newton_iteration;
			
			/* print info about KSP */
			if(debug_print_ksp){
				TRYCXX( KSPView(ksp,PETSC_VIEWER_STDOUT_SELF) );
			}
			
			/* destroy KSP solver */
			TRYCXX( KSPDestroy(&ksp) ); 

			/* compute norm of inner error */
			if(debug_print_it){
				/* compute norm(H*delta-g) */
				TRYCXX( VecDuplicate(g_Vec,&g_inner_Vec) );
				TRYCXX( MatMult(H_petsc,delta_Vec,g_inner_Vec) );
				TRYCXX( VecAXPY(g_inner_Vec,-1.0,g_Vec) );
				
				TRYCXX( VecNorm(g_inner_Vec, NORM_2, &gnorm_inner) );
			}
			
			/* use solution from KSP to update newton iterations */
			this->timer_update.start();
			 TRYCXX( VecAXPY(xk_Vec,this->newton_coeff,delta_Vec)); /* x = x + delta; */
			this->timer_update.stop();
			
			/* recompute integrals, gradient, function value */
			this->timer_integrate.start();
			 compute_integrals(integralsk_Vec, xk_Vec, true);
			this->timer_integrate.stop();
			this->timer_g.start();
			 compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_g.stop();
			this->timer_fs.start();
			 fx = compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_fs.stop();			
			
			it++;

			/* print progress of algorithm */
			if(debug_print_it){
				TRYCXX( VecNorm(delta_Vec, NORM_2, &deltanorm) );
				
				coutMaster << "\033[33m   it = \033[0m" << it;
				std::streamsize ss = std::cout.precision();
				coutMaster << ", \t\033[36mfx = \033[0m" << std::setprecision(17) << fx << std::setprecision(ss);
				coutMaster << ", \t\033[36mnorm(delta) = \033[0m" << std::setprecision(17) << deltanorm << std::setprecision(ss);
				coutMaster << ", \t\033[36mnorm(g_inner) = \033[0m" << gnorm_inner;
				coutMaster << ", \t\033[36mnorm(g_outer) = \033[0m" << gnorm << std::endl;

				/* log function value */
				LOG_FX(fx)
			}


			///* monitor - export values of stopping criteria */
			//if(this->monitor && GlobalManager.get_rank() == 0){
				////TODO: this could be done in a different way
				//std::ofstream myfile;
				//myfile.open("log/entropysolvernewton_monitor.m", std::fstream::in | std::fstream::out | std::fstream::app);

				//std::streamsize ss = myfile.precision();
				//myfile << std::setprecision(17);
			
				//myfile << "fx(" << it << ") = " << fx << "; "; 
				//myfile << "alpha_bb(" << it << ") = " << alpha_bb << "; ";
				//myfile << "gnorm(" << it << ") = " << gnorm << "; ";
				//myfile << std::endl;
			
				//myfile << std::setprecision(ss);
				//myfile.close();			
			//}

		}

		/* store number of iterations */
		this->it_lasts[k] = it;
		this->it_sums[k] += it;
		this->itksp_lasts[k] = itksp;
		this->itksp_sums[k] += itksp;
		this->fxs[k] = fx;
		this->gnorms[k] = gnorm;

		/* -------------- end of SPG algorithm ------------ */

		/* restore subvectors */
		TRYCXX( VecRestoreSubVector(x_Vec, k_is, &xk_Vec) );
		TRYCXX( VecRestoreSubVector(moments_Vec, k_is, &momentsk_Vec) );
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* destroy index set */
		TRYCXX( ISDestroy(&k_is) );
		TRYCXX( ISDestroy(&integralsk_is) );	

	} /* endfor through clusters */


	this->timer_solve.stop(); 

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::compute_moments_data() {
	LOG_FUNC_BEGIN

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec x_gammak_Vec = x_gammak->get_vector();
	Vec x_gammak_power_Vec = x_gammak_power->get_vector();

	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments_data->get_vector();
	
	Vec gammak_Vec;
	IS gammak_is;

	double *moments_arr, mysum, gammaksum;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
	for(int k=0;k<entropydata->get_K();k++){
		/* get gammak */
		this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		/* compute gammaksum */
		TRYCXX( VecSum(gammak_Vec, &gammaksum) );

		TRYCXX( VecPointwiseMult(x_gammak_Vec, gammak_Vec, x_Vec) ); /* (gammak.*x)^1 */
		TRYCXX( VecCopy(x_gammak_Vec, x_gammak_power_Vec) );

		for(int km=0; km < entropydata->get_Km(); km++){
		
			TRYCXX( VecSum(x_gammak_power_Vec, &mysum) );

			/* store computed moment */
			if(gammaksum != 0){
				moments_arr[k*this->entropydata->get_Km() + km] = mysum/gammaksum;
			} else {
				moments_arr[k*this->entropydata->get_Km() + km] = 0.0;
			}
	
			/* compute x_power_gammak */
			TRYCXX( VecPointwiseMult(x_gammak_power_Vec, x_gammak_power_Vec, x_gammak_Vec) ); 
		}

		TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	
	}

	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	LOG_FUNC_END
}


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

template<>
void EntropySolverNewton<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
	LOG_FUNC_BEGIN

	int T = entropydata->get_T();
	int Tlocal = entropydata->get_decomposition()->get_Tlocal();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* lambda vector for Dlib integration */
	column_vector lambda(Km);
    auto mom_function = [&](double x)->double { return gg(x, 0, lambda);};
    double F_;

	/* update gamma_solver data - prepare new linear term */
	/* theta includes all moments */
	const double *lambda_arr;
	TRYCXX( VecGetArrayRead(entropydata->get_lambda()->get_vector(), &lambda_arr) );
	
	const double *x_arr;
	TRYCXX( VecGetArrayRead(entropydata->get_x()->get_vector(), &x_arr) );
	
	double *residuum_arr;
	TRYCXX( VecGetArray(residuum->get_vector(), &residuum_arr) );

	double mysum, x_power;
	for(int k=0;k<K;k++){
		/* compute int_X exp(-sum lambda_j x^j) dx for this cluster
		/* from arr to Dlib-vec */
		for(int km=0;km<Km;km++){
			lambda(km) = lambda_arr[k*Km+km];
		}
		F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
		F_ = log(F_);

		for(int t=0;t<Tlocal;t++){
			/* compute sum_{j=1}^{Km} lambda_j*x^j */
			mysum = 0.0;
			x_power = x_arr[t]; /* x^1 */
			for(int km=0;km<Km;km++){
				mysum += lambda_arr[k*Km+km]*x_power;
				x_power *= x_arr[t]; /* x_power = x^(km+1) */
			}

			residuum_arr[t*K + k] = mysum + F_;
		}
	}

	/* coeffs of A_shared are updated via computation of Theta :) */

	/* restore arrays */
	TRYCXX( VecRestoreArray(residuum->get_vector(), &residuum_arr) );
	TRYCXX( VecRestoreArrayRead(entropydata->get_x()->get_vector(), &x_arr) );
	TRYCXX( VecRestoreArrayRead(entropydata->get_lambda()->get_vector(), &lambda_arr) );
		
	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN

	int Km = entropydata->get_Km();
    
    double *g_arr;
    double *integrals_arr;
    double *moments_arr;

	TRYCXX( VecGetArray(g_Vec, &g_arr) );
	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
    
    /* compute gradient */
    for (int km = 0; km < Km; km++){
        g_arr[km] = moments_arr[km] - integrals_arr[km+1]/integrals_arr[0];
	} 

	TRYCXX( VecRestoreArray(g_Vec, &g_arr) );
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );
	
	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::compute_hessian(Vec &integrals_Vec) {
	LOG_FUNC_BEGIN

	int Km = entropydata->get_Km();
    
    double *integrals_arr;

	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );

    /* fill Hessian matrix */
	for(int km1=0; km1 < Km; km1++){
		for(int km2=0; km2 < Km; km2++){
			double value =  integrals_arr[km1+km2+2]/integrals_arr[0] - integrals_arr[km1+1]*integrals_arr[km2+1]/(integrals_arr[0]*integrals_arr[0]);
			TRYCXX( MatSetValue(H_petsc, km1, km2, value, INSERT_VALUES) );
		}
	} 

	/* assemble matrix */
	TRYCXX( MatAssemblyBegin(H_petsc,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(H_petsc,MAT_FINAL_ASSEMBLY) );	

	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
	
	LOG_FUNC_END
}

template<>
double EntropySolverNewton<PetscVector>::compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN
	
	double f, momTlambda;
    double *integrals_arr;

	TRYCXX( VecDot(moments_Vec, lambda_Vec, &momTlambda) );

	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );
	f = log(integrals_arr[0]) + momTlambda;
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
		
	LOG_FUNC_END

	return f;
}

template<>
void EntropySolverNewton<PetscVector>::compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all) {
	LOG_FUNC_BEGIN

	int Km = entropydata->get_Km();
	int Km_int; /* number of computed integrals */
	if(compute_all){
		Km_int = 2*Km;
	} else {
		Km_int = 0;
	}

	/* create instance of integration tool (if not been created before) */
	if(!this->entropyintegration){
		this->entropyintegration = new EntropyIntegration(Km, Km_int);
	}

	double *integrals_arr;
	double *lambda_arr;
	TRYCXX( VecGetArray(lambda_Vec, &lambda_arr) );
	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );

	//// TODO: temp
    column_vector lambda_Dlib(Km);
    for(int km=0;km<Km;km++){
		lambda_Dlib(km) = lambda_arr[km];
	}

    /* compute integrals */
    for(int km = 0; km<=Km_int;km++){
		auto mom_function = [&](double x)->double { return gg(x, km, lambda_Dlib);};
		integrals_arr[km] = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
	}

	TRYCXX( VecRestoreArray(lambda_Vec, &lambda_arr) );
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );

	LOG_FUNC_END
}


}
} /* end namespace */

#endif
