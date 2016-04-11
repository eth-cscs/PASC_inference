#ifndef PASC_CGQPSOLVER_H
#define	PASC_CGQPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define CGQPSOLVER_DEFAULT_MAXIT 1000;
#define CGQPSOLVER_DEFAULT_EPS 0.0001;
#define CGQPSOLVER_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class CGQPSolverSetting : public QPSolverSetting {
	public:
		CGQPSolverSetting() {
			this->maxit = CGQPSOLVER_DEFAULT_MAXIT;
			this->eps = CGQPSOLVER_DEFAULT_EPS;
			this->debug_mode = CGQPSOLVER_DEFAULT_DEBUG_MODE;
		};
		~CGQPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - maxit:      " << this->maxit << std::endl;
			output <<  " - eps:        " << this->eps << std::endl;
			output <<  " - debug_mode: " << this->debug_mode << std::endl;

		};

		std::string get_name() const {
			return "CG SolverSetting";
		};
		
};


/* CGQPSolver */ 
template<class VectorBase>
class CGQPSolver: public QPSolver<VectorBase> {
	protected:
		const QPData<VectorBase> *qpdata; /* data on which the solver operates */
	
		/* temporary vectors used during the solution process */
		void allocate_temp_vectors();
		void free_temp_vectors();
		GeneralVector<VectorBase> *g; /* gradient */
		GeneralVector<VectorBase> *p; /* A-conjugate vector */
		GeneralVector<VectorBase> *Ap; /* A*p */
	
	public:
		CGQPSolverSetting setting;

		CGQPSolver();
		CGQPSolver(const QPData<VectorBase> &new_qpdata); 
		~CGQPSolver();


		void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		void print(std::ostream &output) const;
		void printcontent(std::ostream &output) const;
		void printstatus(std::ostream &output) const;
		void printtimer(std::ostream &output) const;
		
		std::string get_name() const;



};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(){
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	/* temp vectors */
	g = NULL;
	p = NULL;
	Ap = NULL;
	
	/* iterations counters */
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;	
	
	/* timers */
	
	
	/* function value */
	this->fx = std::numeric_limits<double>::max();
}

template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(const QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();

	/* iterations counters */
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;	
	
	/* timers */


	this->fx = std::numeric_limits<double>::max();
}


/* destructor */
template<class VectorBase>
CGQPSolver<VectorBase>::~CGQPSolver(){
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)DESTRUCTOR" << std::endl;

	/* free temp vectors */
	free_temp_vectors();
}

/* prepare temp_vectors */
template<class VectorBase>
void CGQPSolver<VectorBase>::allocate_temp_vectors(){
	GeneralVector<VectorBase> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<VectorBase>(*pattern);
	p = new GeneralVector<VectorBase>(*pattern);
	Ap = new GeneralVector<VectorBase>(*pattern);	
	
}

/* destroy temp_vectors */
template<class VectorBase>
void CGQPSolver<VectorBase>::free_temp_vectors(){
	free(g);
	free(p);
	free(Ap);
	
}


/* print info about problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)FUNCTION: print" << std::endl;

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	coutMaster.push();
	setting.print(output);
	coutMaster.pop();

	/* print settings */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}
		
}

/* print content of solver */
template<class VectorBase>
void CGQPSolver<VectorBase>::printcontent(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)FUNCTION: printcontent" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
}

template<class VectorBase>
void CGQPSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)FUNCTION: printstatus" << std::endl;

	output <<  this->get_name() << std::endl;
	output <<  " - it:          " << this->it_last << std::endl;
	output <<  " - hess mult:   " << this->hessmult_last << std::endl;
	output <<  " - fx:          " << this->fx << std::endl;	
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

}

template<class VectorBase>
void CGQPSolver<VectorBase>::printtimer(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
/*	output <<  " - it all =       " << this->it_sum << std::endl;
	output <<  " - hessmult all = " << this->hessmult_sum << std::endl;
	output <<  " - timers all" << std::endl;
	output <<  "  - t_solve =      " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_project =    " << this->timer_projection.get_value_sum() << std::endl;
	output <<  "  - t_matmult =    " << this->timer_matmult.get_value_sum() << std::endl;
	output <<  "  - t_dot =        " << this->timer_dot.get_value_sum() << std::endl;
	output <<  "  - t_update =     " << this->timer_update.get_value_sum() << std::endl;
	output <<  "  - t_stepsize =   " << this->timer_stepsize.get_value_sum() << std::endl;
	output <<  "  - t_fs =         " << this->timer_fs.get_value_sum() << std::endl;
	output <<  "  - t_other =      " << this->timer_solve.get_value_sum() - (this->timer_projection.get_value_sum() + this->timer_matmult.get_value_sum() + this->timer_dot.get_value_sum() + this->timer_update.get_value_sum() + this->timer_stepsize.get_value_sum() + this->timer_fs.get_value_sum()) << std::endl;
*/
}


template<class VectorBase>
std::string CGQPSolver<VectorBase>::get_name() const {
	return "Conjugate Gradient method for QP";
}


/* solve the problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)FUNCTION: solve" << std::endl;

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
	pVector p = *(this->p); /* A-conjugate vector */
	pVector Ap = *(this->Ap); /* A*p */

	x = x0; /* set approximation as initial */

	int it = 0; /* iteration counter */
	int hessmult = 0; /* number of hessian multiplications */
	double normg, alpha, beta, pAp, gg, gg_old;
	
	g = A*x; hessmult += 1; g -= b; /* compute gradient */
	p = g; /* initial conjugate direction */

	gg = dot(g,g);
	normg = std::sqrt(gg);

	while(normg > this->setting.eps && it < this->setting.maxit){
		/* compute new approximation */

		Ap = A*p; hessmult += 1;

		/* compute step-size */			
		pAp = dot(Ap,p);
		alpha = gg/pAp;

		/* set new approximation, x = x - alpha*p */
		x -= alpha*p; 

		/* compute gradient recursively, g = g - alpha*Ap */
		g -= alpha*Ap; 
		gg_old = gg;
		gg = dot(g,g);
		normg = std::sqrt(gg);
			
		/* compute new A-orthogonal vector, p = g + beta*p */
		beta = gg/gg_old;
		p *= beta;
		p += g;
		
		if(this->setting.debug_mode >= 10){
			coutMaster << "it " << it << ": ||g|| = " << normg << std::endl;
		}

		if(this->setting.debug_mode >= 100){
			coutMaster << "x = " << x << std::endl;
			coutMaster << "g = " << g << std::endl;
			coutMaster << "p = " << p << std::endl;
			coutMaster << "Ap = " << Ap << std::endl;
			coutMaster << "pAp = " << pAp << std::endl;
			coutMaster << "alpha = " << alpha << std::endl;
			coutMaster << "beta = " << beta << std::endl;
			coutMaster << "gg = " << gg << std::endl;
			coutMaster << "gg_old = " << gg_old << std::endl;
			coutMaster << "normg = " << normg << std::endl;

			coutMaster << "------------------------------------" << std::endl;
		}

				
		it += 1;

	}
		
	/* print output */
	if(this->setting.debug_mode >= 10){
		coutMaster << "------------------------" << std::endl;
		coutMaster << " it_cg = " << it << std::endl;
		coutMaster << " norm_g = " << normg << std::endl;
		coutMaster << " hessmult = " << hessmult << std::endl;
	}

	this->it_sum += it;
	this->hessmult_sum += hessmult;
	this->it_last = it;
	this->hessmult_last = hessmult;

	this->fx = normg; /* fx = norm(g) */
	
}

template<class VectorBase>
double CGQPSolver<VectorBase>::get_fx() const {
	if(setting.debug_mode >= 11) coutMaster << "(CGQPSolver)FUNCTION: get_fx()" << std::endl;
	
	return this->fx;	
}

template<class VectorBase>
int CGQPSolver<VectorBase>::get_it() const {
	return this->it_last;
}

template<class VectorBase>
int CGQPSolver<VectorBase>::get_hessmult() const {
	return this->hessmult_last;
}



} /* end namespace */

#endif
