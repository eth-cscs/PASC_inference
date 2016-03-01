#ifndef CGQPSOLVER_H
#define	CGQPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "data/qpdata.h"
#include "result/qpresult.h"

#define CGQPSOLVER_DEFAULT_MAXIT 1000;
#define CGQPSOLVER_DEFAULT_EPS 0.0001;

namespace pascinference {

/* settings */
class CGQPSolverSetting : public QPSolverSetting {
	protected:
		int maxit;
		double eps;
		
	public:
		CGQPSolverSetting() {
			maxit = CGQPSOLVER_DEFAULT_MAXIT;
			eps = CGQPSOLVER_DEFAULT_EPS;
		};
		~CGQPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  CGQPSolverSettings:" << std::endl;
			output << "   - maxit: " << maxit << std::endl;
			output << "   - eps: " << eps << std::endl;

		};
		
};


/* CGQPSolver */ 
template<class VectorBase>
class CGQPSolver: public QPSolver<VectorBase> {
	protected:
		const QPData<VectorBase> *data; /* data on which the solver operates */
		const QPResult<VectorBase> *result; /* here solver stores results */
	
		/* temporary vectors used during the solution process */
		void allocate_temp_vectors();
		void free_temp_vectors();
		GeneralVector<VectorBase> *g; /* gradient */
		GeneralVector<VectorBase> *p; /* A-conjugate vector */
		GeneralVector<VectorBase> *Ap; /* A*p */
	
	public:
		CGQPSolverSetting setting;

		CGQPSolver();
		CGQPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result); 
		~CGQPSolver();


		void solve();
		void solve(SolverType type){};

		void print(std::ostream &output) const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(CGQPSolver)CONSTRUCTOR" << std::endl;

	data = NULL;
	result = NULL;
	
	/* temp vectors */
	g = NULL;
	p = NULL;
	Ap = NULL;
}

template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result){
	data = &new_data;
	result = &new_result;
	
	/* allocate temp vectors */
	allocate_temp_vectors();
	
}


/* destructor */
template<class VectorBase>
CGQPSolver<VectorBase>::~CGQPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(CGQPSolver)DESTRUCTOR" << std::endl;

	/* free temp vectors */
	free_temp_vectors();
}

/* prepare temp_vectors */
template<class VectorBase>
void CGQPSolver<VectorBase>::allocate_temp_vectors(){
	GeneralVector<VectorBase> *pattern = data->b; /* I will allocate temp vectors subject to linear term */

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
	if(DEBUG_MODE >= 100) std::cout << "(CGQPSolver)FUNCTION: print" << std::endl;

	output << " CGQPSolver" << std::endl;
	
	/* print settings */
	output << setting;
		
}

/* solve the problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::solve() {
	if(DEBUG_MODE >= 100) std::cout << "(CGQPSolver)FUNCTION: solve" << std::endl;

	/* pointers to data */
	GeneralMatrix<VectorBase> *A = data->A;
	GeneralVector<VectorBase> *b = data->b;
	GeneralVector<VectorBase> *x0 = data->x0;

	/* pointers to result */
	GeneralVector<VectorBase> *x = result->x;

	(*x) = (*x0); /* set approximation as initial */

	/* CG method */
	GeneralVector<VectorBase> g(*b); /* gradient */
	GeneralVector<VectorBase> p(*b); /* A-conjugate vector */
	GeneralVector<VectorBase> Ap(*b); /* A*p */

	int it = 0; /* iteration counter */
	int hess_mult = 0; /* number of hessian multiplications */
	double normg, alpha, beta, pAp, gg, gg_old;
	
	g = (*A)*(*x); hess_mult += 1; g -= (*b); /* compute gradient */
	p = g; /* initial conjugate gradient */

	gg = dot(g,g);
	normg = std::sqrt(gg);

	while(normg > 0.001 && it < 10000){
		/* compute new approximation */

		Ap = (*A)*p; hess_mult += 1;
			
		pAp = dot(Ap,p);
		alpha = gg/pAp; /* compute step-size */
		(*x) -= alpha*p; /* set new approximation */

		/* compute gradient recursively */
		g -= alpha*Ap; 
		gg_old = gg;
		gg = dot(g,g);
		normg = std::sqrt(gg);
			
		/* compute new A-orthogonal vector */
		beta = gg/gg_old;
		p *= beta;
		p += g;
		
		std::cout << "it " << it << ": ||g|| = " << normg << std::endl;

		if(DEBUG_MODE >= 10){
			std::cout << "x = " << *x << std::endl;
			std::cout << "g = " << g << std::endl;
			std::cout << "p = " << p << std::endl;
			std::cout << "Ap = " << Ap << std::endl;
			std::cout << "pAp = " << pAp << std::endl;
			std::cout << "alpha = " << alpha << std::endl;
			std::cout << "beta = " << beta << std::endl;
			std::cout << "gg = " << gg << std::endl;
			std::cout << "gg_old = " << gg_old << std::endl;
			std::cout << "normg = " << normg << std::endl;

			std::cout << "------------------------------------" << std::endl;
		}

				
		it += 1;

	}
		
	/* print output */
	std::cout << "------------------------" << std::endl;
	std::cout << " it_cg = " << it << std::endl;
	std::cout << " norm_g = " << normg << std::endl;
	std::cout << " hess_mult = " << hess_mult << std::endl;


	
}


} /* end namespace */

#endif
