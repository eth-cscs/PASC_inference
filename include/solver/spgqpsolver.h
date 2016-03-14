#ifndef PASC_SPGQPSOLVER_H
#define	PASC_SPGQPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include <list>
#include <algorithm>

#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define SPGQPSOLVER_DEFAULT_MAXIT 1000;
#define SPGQPSOLVER_DEFAULT_EPS 0.0001;

#define SPGQPSOLVER_DEFAULT_M 10;
#define SPGQPSOLVER_DEFAULT_GAMMA 0.9;
#define SPGQPSOLVER_DEFAULT_SIGMA2 1.0;
#define SPGQPSOLVER_DEFAULT_ALPHAINIT 2.0;

namespace pascinference {

/* settings */
class SPGQPSolverSetting : public QPSolverSetting {
	public:
		int m; /* size of fs */
		double gamma; 
		double sigma2;
		double alphainit; /* initial step-size */

		SPGQPSolverSetting() {
			maxit = SPGQPSOLVER_DEFAULT_MAXIT;
			eps = SPGQPSOLVER_DEFAULT_EPS;
			debug_mode = DEBUG_MODE;

			m = SPGQPSOLVER_DEFAULT_M;
			gamma = SPGQPSOLVER_DEFAULT_GAMMA;
			sigma2 = SPGQPSOLVER_DEFAULT_SIGMA2;
			alphainit = SPGQPSOLVER_DEFAULT_ALPHAINIT;

		};
		~SPGQPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  SPGQPSolverSettings:" << std::endl;
			output << "   - maxit:      " << maxit << std::endl;
			output << "   - eps:        " << eps << std::endl;
			output << "   - debug_mode: " << debug_mode << std::endl;

			output << "   - m:          " << m << std::endl;
			output << "   - gamma:      " << gamma << std::endl;
			output << "   - sigma2:     " << sigma2 << std::endl;
			output << "   - alphainit:  " << alphainit << std::endl;

		};
		
};


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
		
		friend std::ostream &operator<<(std::ostream &output, SPGQPSolver_fs fs);
};



/* SPGQPSolver */ 
template<class VectorBase>
class SPGQPSolver: public QPSolver<VectorBase> {
	private:
		Timer timer_solve; /* total solution time of SPG algorithm */
		Timer timer_projection; /* the sum of time necessary to perform projections */
		Timer timer_matmult; /* the sum of time necessary to perform matrix multiplication */
		Timer timer_dot; /* the sum of time necessary to compute dot_products */
		Timer timer_update; /* total time of vector updates */
		Timer timer_stepsize; /* total time of step-size computation */
		Timer timer_fs; /* total time of manipulation with fs vector during iterations */

		QPData<VectorBase> *qpdata; /* data on which the solver operates */
	
		/* temporary vectors used during the solution process */
		void allocate_temp_vectors();
		void free_temp_vectors();
		GeneralVector<VectorBase> *g; /* gradient */
		GeneralVector<VectorBase> *d; /* projected gradient */
		GeneralVector<VectorBase> *Ad; /* A*d */
		GeneralVector<VectorBase> *temp; /* general temp vector */

	public:
		SPGQPSolverSetting setting;

		SPGQPSolver();
		SPGQPSolver(QPData<VectorBase> &new_qpdata); 
		~SPGQPSolver();

		void solve();
		void solve(SolverType type){};
		double get_fx() const;
		int get_it() const;

		void print(std::ostream &output) const;
		void printstatus(std::ostream &output) const;		
		void printtimer(std::ostream &output) const;

		std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(){
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	/* temp vectors */
	g = NULL;
	d = NULL;
	Ad = NULL;
	temp = NULL;

	this->it = 0;
	this->hessmult = 0;
	
	this->fx = std::numeric_limits<double>::max();
	
	/* prepare timers */
	this->timer_projection.restart();
	this->timer_matmult.restart();
	this->timer_dot.restart();
	this->timer_update.restart();
	this->timer_stepsize.restart();
	this->timer_fs.restart();
	this->timer_solve.restart();	
	
}

template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();

	this->it = 0;
	this->hessmult = 0;

	this->fx = std::numeric_limits<double>::max();

	/* prepare timers */
	this->timer_projection.restart();
	this->timer_matmult.restart();
	this->timer_dot.restart();
	this->timer_update.restart();
	this->timer_stepsize.restart();
	this->timer_fs.restart();
	this->timer_solve.restart();	

}


/* destructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::~SPGQPSolver(){
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)DESTRUCTOR" << std::endl;

	/* free temp vectors */
	free_temp_vectors();
}

/* prepare temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::allocate_temp_vectors(){
	GeneralVector<VectorBase> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<VectorBase>(*pattern);
	d = new GeneralVector<VectorBase>(*pattern);
	Ad = new GeneralVector<VectorBase>(*pattern);	
	temp = new GeneralVector<VectorBase>(*pattern);	
	
}

/* destroy temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::free_temp_vectors(){
	free(g);
	free(d);
	free(Ad);
	free(temp);
	
}


/* print info about problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)FUNCTION: print" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;
		
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)FUNCTION: printstatus" << std::endl;

	output << this->get_name() << std::endl;
	output << " - it:          " << this->it << std::endl;
	output << " - fx:          " << this->fx << std::endl;	
	output << " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;
}

template<class VectorBase>
void SPGQPSolver<VectorBase>::printtimer(std::ostream &output) const {
	output << this->get_name() << std::endl;
	output << "    - it =        " << this->it << std::endl;
	output << "    - hessmult =  " << this->hessmult << std::endl;
	output << "    - timers" << std::endl;
	output << "     - t_solve =  " << this->timer_solve.get_value_sum() << std::endl;
	output << "     - t_project =  " << this->timer_projection.get_value_sum() << std::endl;
	output << "     - t_matmult =  " << this->timer_matmult.get_value_sum() << std::endl;
	output << "     - t_dot =      " << this->timer_dot.get_value_sum() << std::endl;
	output << "     - t_update =   " << this->timer_update.get_value_sum() << std::endl;
	output << "     - t_stepsize = " << this->timer_stepsize.get_value_sum() << std::endl;
	output << "     - t_fs =       " << this->timer_fs.get_value_sum() << std::endl;
//	output << "     - t_other =    ", this->timer_solve.get_value_sum() - (this->timer_projection.get_value_sum() + this->timer_matmult.get_value_sum() + this->timer_dot.get_value_sum() + this->timer_update.get_value_sum() + this->timer_stepsize.get_value_sum() + this->timer_fs.get_value_sum()) << std::endl;
}


template<class VectorBase>
std::string SPGQPSolver<VectorBase>::get_name() const {
	return "Spectral Projected Gradient method for QP";
}

/* solve the problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::solve() {
	if(setting.debug_mode >= 100) std::cout << "(SPGQPSolver)FUNCTION: solve" << std::endl;

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
	SPGQPSolver_fs fs(setting.m); /* store function values for generalized Armijo condition */
	double fx_max; /* max(fs) */
	double xi, beta_bar, beta_hat, beta; /* for Armijo condition */
	double dd; /* dot(d,d) */
	double gd; /* dot(g,d) */
	double dAd; /* dot(Ad,d) */
	double alpha_bb; /* BB step-size */

	/* initial step-size */
	alpha_bb = setting.alphainit;

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
	 this->fx = fx;
	 fs.init(fx);	
	this->timer_fs.stop();

	/* main cycle */
	while(it < setting.maxit){

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

		/* dd = dot(d,d) */
		/* dAd = dot(Ad,d) */
		/* gd = dot(g,d) */
		this->timer_dot.start();
		 dd = dot(d,d);
		 dAd = dot(Ad,d);
		 gd = dot(g,d);
		this->timer_dot.stop();

		/* stopping criteria */
		if(dd < setting.eps){
			break;
		}

		/* fx_max = max(fs) */
		this->timer_fs.start();
		 fx_max = fs.get_max();	
		this->timer_fs.stop();
		
		/* compute step-size from A-condition */
		this->timer_stepsize.start();
		 xi = (fx_max - fx)/dAd;
		 beta_bar = -gd/dAd;
		 beta_hat = setting.gamma*beta_bar + std::sqrt(setting.gamma*setting.gamma*beta_bar*beta_bar + 2*xi);

		 /* beta = min(sigma2,beta_hat) */
		 if(beta_hat < setting.sigma2){
			beta = beta_hat;
		 } else {
			beta = setting.sigma2;
		 }
		this->timer_stepsize.stop();

		/* update approximation and gradient */
		this->timer_update.start();
		 x += beta*d; /* x = x + beta*d */
		 g += beta*Ad; /* g = g + beta*Ad */
		this->timer_update.stop();

		/* compute new function value using gradient and update fs list */
		this->timer_fs.start();
		 fx = get_fx();
		 fs.update(fx);
		this->timer_fs.stop();

		/* update BB step-size */
		this->timer_stepsize.start();
		 alpha_bb = dd/dAd;
		this->timer_stepsize.stop();

		/* print qpdata */
		if(setting.debug_mode >= 10){
			std::cout << "x: " << x << std::endl;
			std::cout << "d: " << d << std::endl;
			std::cout << "g: " << g << std::endl;
			std::cout << "Ad: " << Ad << std::endl;
			
		}

		/* print progress of algorithm */
		if(setting.debug_mode >= 3){
			std::cout << "\033[33m   it = \033[0m" << it;
			std::cout << ", \t\033[36mfx = \033[0m" << fx;
			std::cout << ", \t\033[36mdd = \033[0m" << dd << std::endl;
		}

		if(setting.debug_mode >= 5){
			std::cout << "\033[36m    alpha_bb = \033[0m" << alpha_bb << ",";
			std::cout << "\033[36m dAd = \033[0m" << dAd << ",";
			std::cout << "\033[36m gd = \033[0m" << gd << std::endl;
			
			std::cout << "\033[36m    fx = \033[0m" << fx << ",";
			std::cout << "\033[36m fx_max = \033[0m" << fx_max << ",";
			std::cout << "\033[36m xi = \033[0m" << xi << std::endl;
			
			std::cout << "\033[36m    beta_bar = \033[0m" << beta_bar << ",";
			std::cout << "\033[36m beta_hat = \033[0m" << beta_hat << ",";
			std::cout << "\033[36m beta = \033[0m" << beta << std::endl;
			
		}
		
		/* increase iteration counter */
		it += 1;
	} /* main cycle end */

	this->it += it;
	this->hessmult += hessmult;
	this->fx = fx;
	this->timer_solve.stop();

	/* very short info */
	if(setting.debug_mode >= 3){
		Message_info_value("   - it    = ",it);
		Message_info_time("   - time  = ",this->timer_solve.get_value_last());

	}

	
}

/* compute function value using inner *x and already computed *g */
template<class VectorBase>
double SPGQPSolver<VectorBase>::get_fx() const {
	if(setting.debug_mode >= 11) std::cout << "(SPGQPSolver)FUNCTION: get_fx()" << std::endl;
	
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
	fx = 0.5*dot(temp,x);

	return fx;	
}

template<class VectorBase>
int SPGQPSolver<VectorBase>::get_it() const {
	return this->it;
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
std::ostream &operator<<(std::ostream &output, SPGQPSolver_fs fs)
{
	int j, list_size;
	std::list<double>::iterator it; /* iterator through list */

	it = fs.fs_list.begin();
	list_size = fs.fs_list.size(); /* = m? */
	
	output << "[ ";
	/* for each component go throught the list */
	for(j=0;j<list_size;j++){
		output << *it;
		if(j < list_size-1){ 
				/* this is not the last node */
				output << ", ";
				it++;
		}
	}
	output << " ]";
			
	return output;
}



} /* end namespace */

#endif