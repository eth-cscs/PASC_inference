#ifndef PASC_CGQPSOLVER_H
#define	PASC_CGQPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define CGQPSOLVER_DEFAULT_MAXIT 1000;
#define CGQPSOLVER_DEFAULT_EPS 0.0001;

namespace pascinference {

/* settings */
class CGQPSolverSetting : public QPSolverSetting {
	public:
		CGQPSolverSetting() {
			this->maxit = CGQPSOLVER_DEFAULT_MAXIT;
			this->eps = CGQPSOLVER_DEFAULT_EPS;
			this->debug_mode = DEBUG_MODE;
		};
		~CGQPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << offset << this->get_name() << std::endl;
			output << offset << " - maxit:      " << this->maxit << std::endl;
			output << offset << " - eps:        " << this->eps << std::endl;
			output << offset << " - debug_mode: " << this->debug_mode << std::endl;

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
		void solve(SolverType type){};

		void print(std::ostream &output) const;
		std::string get_name() const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(){
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(CGQPSolver)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	/* temp vectors */
	g = NULL;
	p = NULL;
	Ap = NULL;
	
	this->fx = std::numeric_limits<double>::max();
}

template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(const QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();

	this->fx = std::numeric_limits<double>::max();
}


/* destructor */
template<class VectorBase>
CGQPSolver<VectorBase>::~CGQPSolver(){
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(CGQPSolver)DESTRUCTOR" << std::endl;

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
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(CGQPSolver)FUNCTION: print" << std::endl;

	output << offset << this->get_name() << std::endl;
	
	/* print settings */
	offset.push();
	setting.print(output);
	offset.pop();

	/* print settings */
	if(qpdata){
		offset.push();
		qpdata->print(output);
		offset.pop();
	}
		
}

template<class VectorBase>
std::string CGQPSolver<VectorBase>::get_name() const {
	return "Conjugate Gradient method for QP";
}


/* solve the problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::solve() {
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(CGQPSolver)FUNCTION: solve" << std::endl;

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
	int hess_mult = 0; /* number of hessian multiplications */
	double normg, alpha, beta, pAp, gg, gg_old;
	
	g = A*x; hess_mult += 1; g -= b; /* compute gradient */
	p = g; /* initial conjugate direction */

	gg = dot(g,g);
	normg = std::sqrt(gg);

	while(normg > this->setting.eps && it < this->setting.maxit){
		/* compute new approximation */

		Ap = A*p; hess_mult += 1;

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
			coutMaster << offset <<"it " << it << ": ||g|| = " << normg << std::endl;
		}

		if(this->setting.debug_mode >= 100){
			coutMaster << offset <<"x = " << x << std::endl;
			coutMaster << offset <<"g = " << g << std::endl;
			coutMaster << offset <<"p = " << p << std::endl;
			coutMaster << offset <<"Ap = " << Ap << std::endl;
			coutMaster << offset <<"pAp = " << pAp << std::endl;
			coutMaster << offset <<"alpha = " << alpha << std::endl;
			coutMaster << offset <<"beta = " << beta << std::endl;
			coutMaster << offset <<"gg = " << gg << std::endl;
			coutMaster << offset <<"gg_old = " << gg_old << std::endl;
			coutMaster << offset <<"normg = " << normg << std::endl;

			coutMaster << offset <<"------------------------------------" << std::endl;
		}

				
		it += 1;

	}
		
	/* print output */
	if(this->setting.debug_mode >= 10){
		coutMaster << offset <<"------------------------" << std::endl;
		coutMaster << offset <<" it_cg = " << it << std::endl;
		coutMaster << offset <<" norm_g = " << normg << std::endl;
		coutMaster << offset <<" hess_mult = " << hess_mult << std::endl;
	}

	
}


} /* end namespace */

#endif
