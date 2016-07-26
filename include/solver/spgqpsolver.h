#ifndef PASC_SPGQPSOLVER_H
#define	PASC_SPGQPSOLVER_H

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define SPGQPSOLVER_DEFAULT_MAXIT 500
#define SPGQPSOLVER_DEFAULT_EPS 1e-9
#define SPGQPSOLVER_DEFAULT_DEBUG_MODE 0

#define SPGQPSOLVER_DEFAULT_M 20
#define SPGQPSOLVER_DEFAULT_GAMMA 0.9
#define SPGQPSOLVER_DEFAULT_SIGMA1 0.000
#define SPGQPSOLVER_DEFAULT_SIGMA2 1.0
#define SPGQPSOLVER_DEFAULT_ALPHAINIT 2.0
#define SPGQPSOLVER_DEFAULT_DOTFLOOR 2

#define SPGQPSOLVER_STOP_NORMGP true
#define SPGQPSOLVER_STOP_ANORMGP false
#define SPGQPSOLVER_STOP_NORMGP_NORMB false
#define SPGQPSOLVER_STOP_ANORMGP_NORMB false
#define SPGQPSOLVER_STOP_DIFFF false

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

namespace pascinference {

/* for manipulation with fs - function values for generalized Armijo condition */
class SPGQPSolver_fs {
	private:
		int m; /* the length of list */
		std::list<double> fs_list; /* the list with function values */

	public: 
		SPGQPSolver_fs(int new_m);
		void init(double fx);
		double get_max();		
		int get_size();
		void update(double new_fx);
		
		void print(ConsoleOutput &output);
};



/* SPGQPSolver */ 
template<class VectorBase>
class SPGQPSolver: public QPSolver<VectorBase> {
	private:
		Timer timer_solve; 			/**< total solution time of SPG algorithm */
		Timer timer_projection;		/**< the sum of time necessary to perform projections */
		Timer timer_matmult; 		/**< the sum of time necessary to perform matrix multiplication */
		Timer timer_dot; 			/**< the sum of time necessary to compute dot_products */
		Timer timer_update; 		/**< total time of vector updates */
		Timer timer_stepsize;	 	/**< total time of step-size computation */
		Timer timer_fs; 			/**< total time of manipulation with fs vector during iterations */

		/* settings */
		bool stop_normgp; 			/**< stopping criteria based on norm of gP */
		bool stop_Anormgp; 			/**< stopping criteria based on A-norm of gP */
		bool stop_normgp_normb;		/**< stopping criteria based on norm of gP and norm of b */
		bool stop_Anormgp_normb;	/**< stopping criteria based on A-norm of gP and norm of b */
		bool stop_difff;			/**< stopping criteria based on size of decrease of f */

		int m;						/**< size of fs */
		double gamma; 
		double sigma1;
		double sigma2;
		double alphainit;			/** initial step-size */

		QPData<VectorBase> *qpdata; /**< data on which the solver operates */
		double gP; 					/**< norm of projected gradient */
	
		/* temporary vectors used during the solution process */
		void allocate_temp_vectors();
		void free_temp_vectors();
		GeneralVector<VectorBase> *g; 		/**< gradient */
		GeneralVector<VectorBase> *d; 		/**< projected gradient */
		GeneralVector<VectorBase> *Ad; 		/**< A*d */
		GeneralVector<VectorBase> *temp;	/**< general temp vector */

		void compute_dots(double *dd, double *dAd, double *gd) const;

		/* for manipulation with mdot */
		double *Mdots_val;

		#ifdef USE_PETSCVECTOR
			Vec *Mdots_vec;
		#endif

		void set_settings_from_console();
		
	public:
		SPGQPSolver();
		SPGQPSolver(QPData<VectorBase> &new_qpdata); 
		~SPGQPSolver();

		void solve();

		double get_fx() const;
		double get_fx(double fx_old, double beta, double gd, double dAd) const;
		int get_it() const;
		int get_hessmult() const;

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printtimer(ConsoleOutput &output) const;
		void printshort(std::ostringstream &header, std::ostringstream &values) const;

		std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<class VectorBase>
void SPGQPSolver<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("spgqpsolver_maxit", &this->maxit, SPGQPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("spgqpsolver_eps", &this->eps, SPGQPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("spgqpsolver_debug_mode", &this->debug_mode, SPGQPSOLVER_DEFAULT_DEBUG_MODE);
	
	consoleArg.set_option_value("spgqpsolver_m", &this->m, SPGQPSOLVER_DEFAULT_M);	
	consoleArg.set_option_value("spgqpsolver_gamma", &this->gamma, SPGQPSOLVER_DEFAULT_GAMMA);	
	consoleArg.set_option_value("spgqpsolver_sigma1", &this->sigma1, SPGQPSOLVER_DEFAULT_SIGMA1);	
	consoleArg.set_option_value("spgqpsolver_sigma2", &this->sigma2, SPGQPSOLVER_DEFAULT_SIGMA2);	
	consoleArg.set_option_value("spgqpsolver_alphainit", &this->alphainit, SPGQPSOLVER_DEFAULT_ALPHAINIT);	

	consoleArg.set_option_value("spgqpsolver_stop_normgp", &this->stop_normgp, SPGQPSOLVER_STOP_NORMGP);
	consoleArg.set_option_value("spgqpsolver_stop_Anormgp", &this->stop_Anormgp, SPGQPSOLVER_STOP_ANORMGP);
	consoleArg.set_option_value("spgqpsolver_stop_normgp_normb", &this->stop_normgp_normb, SPGQPSOLVER_STOP_NORMGP_NORMB);
	consoleArg.set_option_value("spgqpsolver_stop_Anormgp_normb", &this->stop_Anormgp_normb, SPGQPSOLVER_STOP_ANORMGP_NORMB);
	consoleArg.set_option_value("spgqpsolver_stop_difff", &this->stop_difff, SPGQPSOLVER_STOP_DIFFF);	

}


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(){
	LOG_FUNC_BEGIN

	qpdata = NULL;
	
	/* temp vectors */
	this->g = NULL;
	this->d = NULL;
	this->Ad = NULL;
	this->temp = NULL;

	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;
	
	this->fx = std::numeric_limits<double>::max();
	this->gP = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_projection.restart();
	this->timer_matmult.restart();
	this->timer_dot.restart();
	this->timer_update.restart();
	this->timer_stepsize.restart();
	this->timer_fs.restart();

	LOG_FUNC_END
}

template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();

	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;

	this->fx = std::numeric_limits<double>::max();
	this->gP = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_projection.restart();
	this->timer_matmult.restart();
	this->timer_dot.restart();
	this->timer_update.restart();
	this->timer_stepsize.restart();
	this->timer_fs.restart();
	this->timer_solve.restart();	

	LOG_FUNC_END
}


/* destructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::~SPGQPSolver(){
	LOG_FUNC_BEGIN

	/* free temp vectors */
	free_temp_vectors();

	LOG_FUNC_END
}

/* prepare temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<VectorBase> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<VectorBase>(*pattern);
	d = new GeneralVector<VectorBase>(*pattern);
	Ad = new GeneralVector<VectorBase>(*pattern);	
	temp = new GeneralVector<VectorBase>(*pattern);	

	//TODO: prepare Mdot without petsc?
	/* prepare for Mdot */
	#ifdef USE_PETSCVECTOR
		/* without cuda */
		TRY( PetscMalloc1(3,&Mdots_val) );
		TRY( PetscMalloc1(3,&Mdots_vec) );

		Mdots_vec[0] = d->get_vector();
		Mdots_vec[1] = Ad->get_vector();
		Mdots_vec[2] = g->get_vector();

	#endif
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(g);
	free(d);
	free(Ad);
	free(temp);
	
	#ifdef USE_PETSCVECTOR
		TRY( PetscFree(Mdots_val) );
		TRY( PetscFree(Mdots_vec) );
	#endif
	
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debug_mode: " << this->debug_mode << std::endl;

	output <<  " - m:          " << m << std::endl;
	output <<  " - gamma:      " << gamma << std::endl;
	output <<  " - sigma1:     " << sigma1 << std::endl;
	output <<  " - sigma2:     " << sigma2 << std::endl;
	output <<  " - alphainit:  " << alphainit << std::endl;
	
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}

	output.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_local <<  " - maxit:      " << this->maxit << std::endl;
	output_local <<  " - eps:        " << this->eps << std::endl;
	output_local <<  " - debug_mode: " << this->debug_mode << std::endl;

	output_local <<  " - m:          " << m << std::endl;
	output_local <<  " - gamma:      " << gamma << std::endl;
	output_local <<  " - sigma1:     " << sigma1 << std::endl;
	output_local <<  " - sigma2:     " << sigma2 << std::endl;
	output_local <<  " - alphainit:  " << alphainit << std::endl;

	output_local.synchronize();
	
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output_global, output_local);
		coutMaster.pop();
	}

	LOG_FUNC_END
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  " - it: " << std::setw(6) << this->it_last << ", ";
	output <<  "hess mult: " << std::setw(6) << this->hessmult_last << ", ";
	output <<  "fx: " << std::setw(10) << this->fx << ", ";	
	output <<  "norm(gP): " << std::setw(10) << this->gP << ", ";
	output <<  "used memory: " << std::setw(6) << MemoryCheck::get_virtual() << "%" << std::endl;

	output << " - ";
//	output <<  "t_solve = " << std::setw(10) << this->timer_solve.get_value_last() << ", ";
	output <<  "t_project = " << std::setw(10) << this->timer_projection.get_value_last() << ", ";
	output <<  "t_matmult = " << std::setw(10) << this->timer_matmult.get_value_last() << ", ";
//	output <<  "t_dot = " << std::setw(10) << this->timer_dot.get_value_last() << ", ";
//	output <<  "t_update = " << std::setw(10) << this->timer_update.get_value_last() << ", ";
//	output <<  "t_stepsize = " << std::setw(10) << this->timer_stepsize.get_value_last() << ", ";
//	output <<  "t_fs = " << std::setw(10) << this->timer_fs.get_value_last() << ", ";
	output <<  "t_other = " << std::setw(10) << this->timer_solve.get_value_last() - (this->timer_projection.get_value_last() + this->timer_matmult.get_value_last() + this->timer_dot.get_value_last() + this->timer_update.get_value_last() + this->timer_stepsize.get_value_last() + this->timer_fs.get_value_last()) << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);

	std::streamsize ss = std::cout.precision();

	output <<  "      - fx:           " << std::setw(25) << std::setprecision(17) << this->fx << std::endl;
	output <<  "      - fx_control:   " << std::setw(25) << fx_quadratic+fx_linear << ", log: " << std::setw(25) << log(fx_quadratic+fx_linear)/log(10) << std::endl;
	output <<  "      - fx_linear:    " << std::setw(25) << fx_linear << ", log: " << std::setw(25) << log(fx_linear)/log(10) << std::endl;
	output <<  "      - fx_quadratic: " << std::setw(25) << fx_quadratic << ", log: " << std::setw(25) << log(fx_quadratic)/log(10) << std::endl;
	output <<  "      - norm(gP):     " << std::setw(25) << this->gP << ", log: " << std::setw(25) << log(this->gP)/log(10) << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it all =        " << this->it_sum << std::endl;
	output <<  " - hessmult all =  " << this->hessmult_sum << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =      " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_project =    " << this->timer_projection.get_value_sum() << std::endl;
	output <<  "  - t_matmult =    " << this->timer_matmult.get_value_sum() << std::endl;
	output <<  "  - t_dot =        " << this->timer_dot.get_value_sum() << std::endl;
	output <<  "  - t_update =     " << this->timer_update.get_value_sum() << std::endl;
	output <<  "  - t_stepsize =   " << this->timer_stepsize.get_value_sum() << std::endl;
	output <<  "  - t_fs =         " << this->timer_fs.get_value_sum() << std::endl;
	output <<  "  - t_other =      " << this->timer_solve.get_value_sum() - (this->timer_projection.get_value_sum() + this->timer_matmult.get_value_sum() + this->timer_dot.get_value_sum() + this->timer_update.get_value_sum() + this->timer_stepsize.get_value_sum() + this->timer_fs.get_value_sum()) << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	header << "SPGQP it, ";
	values << this->it_sum << ", ";

	header << "SPGQP hessmult, ";
	values << this->hessmult_sum << ", ";

	header << "SPGQP t all, ";
	values << this->timer_solve.get_value_sum() << ", ";

	header << "SPGQP t project, ";
	values << this->timer_projection.get_value_sum() << ", ";

	header << "SPGQP t matmult, ";
	values << this->timer_matmult.get_value_sum() << ", ";

	header << "SPGQP t dot, ";
	values << this->timer_dot.get_value_sum() << ", ";

	header << "SPGQP t update, ";
	values << this->timer_update.get_value_sum() << ", ";

	header << "SPGQP t stepsize, ";
	values << this->timer_stepsize.get_value_sum() << ", ";

	header << "SPGQP t fs, ";
	values << this->timer_fs.get_value_sum() << ", ";

	header << "SPGQP t other, ";
	values << this->timer_solve.get_value_sum() - (this->timer_projection.get_value_sum() + this->timer_matmult.get_value_sum() + this->timer_dot.get_value_sum() + this->timer_update.get_value_sum() + this->timer_stepsize.get_value_sum() + this->timer_fs.get_value_sum()) << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
std::string SPGQPSolver<VectorBase>::get_name() const {
	return "SPGQP";
}

/* solve the problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); /* stop this timer in the end of solution */

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());
	pVector x0 = *(qpdata->get_x0());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector g = *(this->g); /* gradient */
	pVector d = *(this->d); /* A-conjugate vector */
	pVector Ad = *(this->Ad); /* A*p */

	int it = 0; /* number of iterations */
	int hessmult = 0; /* number of hessian multiplications */

	double fx; /* function value */
	double fx_old; /* f(x_{it - 1}) */
	SPGQPSolver_fs fs(this->m); /* store function values for generalized Armijo condition */
	double fx_max; /* max(fs) */
	double xi, beta_bar, beta_hat, beta; /* for Armijo condition */
	double dd; /* dot(d,d) */
	double gd; /* dot(g,d) */
	double dAd; /* dot(Ad,d) */
	double alpha_bb; /* BB step-size */
	double normb = norm(b); /* norm of linear term used in stopping criteria */

	/* initial step-size */
	alpha_bb = this->alphainit;

	x = x0; /* set approximation as initial */

	this->timer_projection.start();
	 qpdata->get_feasibleset()->project(x); /* project initial approximation to feasible set */
	this->timer_projection.stop();

	/* compute gradient, g = A*x-b */
	this->timer_matmult.start();
	 g = A*x;
	 hessmult += 1; /* there was muliplication by A */
	this->timer_matmult.stop();
	g -= b;

	/* initialize fs */
	this->timer_fs.start();
	 fx = get_fx();
	 fx_old = std::numeric_limits<double>::max();
	 this->fx = fx;
	 fs.init(fx);
	this->timer_fs.stop();

	/* main cycle */
	while(it < this->maxit){
		/* increase iteration counter */
		it += 1;

		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		this->timer_update.start();
		 d = x - alpha_bb*g;
		this->timer_update.stop();

		/* d = P(d) */
		this->timer_projection.start();
		 qpdata->get_feasibleset()->project(d);
		this->timer_projection.stop();

		/* d = d - x */
		this->timer_update.start();
		 d -= x;
		this->timer_update.stop();

		/* Ad = A*d */
		this->timer_matmult.start();
		 Ad = A*d;
		 hessmult += 1;
		this->timer_matmult.stop();

		this->timer_dot.start();
//		  dd = dot(d,d);
//		  dAd = dot(Ad,d); 
//		  gd = dot(g,d);
		 compute_dots(&dd, &dAd, &gd);
		this->timer_dot.stop();

		/* fx_max = max(fs) */
		this->timer_fs.start();
		 fx_max = fs.get_max();	
		this->timer_fs.stop();
		
		/* compute step-size from A-condition */
		this->timer_stepsize.start();
		 xi = (fx_max - fx)/dAd;
		 beta_bar = -gd/dAd;
		 beta_hat = this->gamma*beta_bar + PetscSqrtReal(this->gamma*this->gamma*beta_bar*beta_bar + 2*xi);

		 /* beta = max(sigma1,min(sigma2,beta_hat)) */
		 if(beta_hat < this->sigma1){
			 beta_hat = this->sigma1;
		 }
		 
		 if(beta_hat < this->sigma2){
			beta = beta_hat;
		 } else {
			beta = this->sigma2;
		 }
		this->timer_stepsize.stop();

		/* update approximation and gradient */
		this->timer_update.start();
		 x += beta*d; /* x = x + beta*d */
		 g += beta*Ad; /* g = g + beta*Ad */
		this->timer_update.stop();

		/* compute new function value using gradient and update fs list */
		this->timer_fs.start();
		 fx_old = fx;
//		 fx = get_fx(fx_old,beta,gd,dAd);
		 fx = get_fx();
		 fs.update(fx);
		this->timer_fs.stop();

		/* update BB step-size */
		this->timer_stepsize.start();
		 alpha_bb = dd/dAd;
		this->timer_stepsize.stop();

		this->gP = dd;
		
		/* print qpdata */
		if(this->debug_mode == 10){
			coutMaster << "x: " << x << std::endl;
			coutMaster << "d: " << d << std::endl;
			coutMaster << "g: " << g << std::endl;
			coutMaster << "Ad: " << Ad << std::endl;
		}

		/* print progress of algorithm */
		if(this->debug_mode == 3){
			coutMaster << "\033[33m   it = \033[0m" << it;
			
			std::streamsize ss = std::cout.precision();
			coutMaster << ", \t\033[36mfx = \033[0m" << std::setprecision(17) << fx << std::setprecision(ss);

			coutMaster << ", \t\033[36mgP = \033[0m" << this->gP;
			coutMaster << ", \t\033[36mdd = \033[0m" << dd << std::endl;
		}

		if(this->debug_mode == 4){
			coutAll << "\033[33m   it = \033[0m" << it << std::endl;
		}


		if(this->debug_mode == 5){
			coutMaster << "\033[36m    alpha_bb = \033[0m" << alpha_bb << ",";
			coutMaster << "\033[36m dAd = \033[0m" << dAd << ",";
			coutMaster << "\033[36m gd = \033[0m" << gd << std::endl;
			
			coutMaster << "\033[36m    fx = \033[0m" << fx << ",";
			coutMaster << "\033[36m fx_max = \033[0m" << fx_max << ",";
			coutMaster << "\033[36m xi = \033[0m" << xi << std::endl;
			
			coutMaster << "\033[36m    beta_bar = \033[0m" << beta_bar << ",";
			coutMaster << "\033[36m beta_hat = \033[0m" << beta_hat << ",";
			coutMaster << "\033[36m beta = \033[0m" << beta << std::endl;
			
		}

		/* stopping criteria */
		if( this->stop_difff && abs(fx - fx_old) < this->eps){
			break;
		}
		if(this->stop_normgp && dd < this->eps){
			break;
		}
		if(this->stop_normgp_normb && dd < this->eps*normb){
			break;
		}
		if(this->stop_Anormgp && dAd < this->eps){
			break;
		}
		if(this->stop_Anormgp_normb && dAd < this->eps*normb){
			break;
		}
		
	} /* main cycle end */

	this->it_sum += it;
	this->hessmult_sum += hessmult;
	this->it_last = it;
	this->hessmult_last = hessmult;

	this->fx = fx;
	this->timer_solve.stop();

	/* write info to log file */
	LOG_IT(it)
	LOG_FX(fx)

	LOG_FUNC_END
}

/* compute function value using inner *x and already computed *g */
template<class VectorBase>
double SPGQPSolver<VectorBase>::get_fx() const {
	LOG_FUNC_BEGIN
	
	double fx = std::numeric_limits<double>::max();

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);

	/* pointers to qpdata */
	pVector g = *(this->g);
	pVector x = *(qpdata->get_x());
	pVector b = *(qpdata->get_b());
	pVector temp = *(this->temp);

	/* use computed gradient in this->g to compute function value */
	temp = g - b;
	double tempt = dot(temp,x);
	fx = 0.5*tempt;

	LOG_FUNC_END
	return fx;	
}

/* compute function value using previously computed values */
template<class VectorBase>
double SPGQPSolver<VectorBase>::get_fx(double fx_old, double beta, double gd, double dAd) const {
	LOG_FUNC_BEGIN
	
	double fx = fx_old + beta*(gd + 0.5*beta*dAd);

	LOG_FUNC_END
	return fx;	
}

/* compute dot products */
template<class VectorBase>
void SPGQPSolver<VectorBase>::compute_dots(double *dd, double *dAd, double *gd) const {
	LOG_FUNC_BEGIN

	*dd = dot(*d,*d);
	*dAd = dot(*Ad,*d);
	*gd = dot(*g,*d);

	LOG_FUNC_END
}

template<class VectorBase>
int SPGQPSolver<VectorBase>::get_it() const {
	return this->it_last;
}

template<class VectorBase>
int SPGQPSolver<VectorBase>::get_hessmult() const {
	return this->hessmult_last;
}


/* ---------- SPGQPSolver_fs -------------- */

/* constructor */
SPGQPSolver_fs::SPGQPSolver_fs(int new_m){
	this->m = new_m;
}

/* init the list with function values using one initial fx */
void SPGQPSolver_fs::init(double fx){
	this->fs_list.resize(this->m, fx);
}

/* get the size of the list */
int SPGQPSolver_fs::get_size(){
	return this->m;
}

/* get the value of max value in the list */
double SPGQPSolver_fs::get_max(){
	std::list<double>::iterator it;
	it = std::max_element(this->fs_list.begin(), this->fs_list.end());
	return *it;
}

/* update the list by new value - pop the first and push the new value (FIFO) */
void SPGQPSolver_fs::update(double new_fx){
	this->fs_list.pop_back();
	this->fs_list.push_front(new_fx);
}

/* print the content of the list */
void SPGQPSolver_fs::print(ConsoleOutput &output)
{
	int j, list_size;
	std::list<double>::iterator it; /* iterator through list */

	it = this->fs_list.begin();
	list_size = this->fs_list.size(); /* = m? */
	
	std::ostringstream temp;
	
	output << "[ ";
	/* for each component go throught the list */
	for(j=0;j<list_size;j++){
		temp << *it;
		output << temp.str();
		temp.str("");
		if(j < list_size-1){ 
				/* this is not the last node */
				output << ", ";
				it++;
		}
	}
	output << " ]";
}


/* ----------- PETSC special stuff ------------- */
#ifdef USE_PETSC
/* compute dot products */
template<>
void SPGQPSolver<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const {
	LOG_FUNC_BEGIN

	TRY(VecMDot( Mdots_vec[0], 3, Mdots_vec, Mdots_val) );

	*dd = Mdots_val[0];
	*dAd = Mdots_val[1];
	*gd = Mdots_val[2];

	LOG_FUNC_END
}
#endif



} /* end namespace */

#endif
