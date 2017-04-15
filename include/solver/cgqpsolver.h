/** @file cgqpsolver.h
 *  @brief Conjugate Gradient method for solving Quadratic Programs
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_CGQPSOLVER_H
#define	PASC_CGQPSOLVER_H

#include "pascinference.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define CGQPSOLVER_DEFAULT_MAXIT 1000
#define CGQPSOLVER_DEFAULT_EPS 0.0001
#define CGQPSOLVER_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class CGQPSolver
 *  \brief Conjugate Gradient method for solving Quadratic Programs
 *
 *  For solving QP without constraints.
*/
template<class VectorBase>
class CGQPSolver: public QPSolver<VectorBase> {
	protected:
		void allocate_temp_vectors(); 
		void free_temp_vectors();

		GeneralVector<VectorBase> *g; /**< auxiliary vector used to store gradient */
		GeneralVector<VectorBase> *p; /**< auxiliary vector used to store A-conjugated vector */
		GeneralVector<VectorBase> *Ap; /**< auxiliary vector used to store A times gradient */
	
	public:
		/** @brief general constructor
		* 
		*/
		CGQPSolver();

		/** @brief constructor based on provided data of problem
		* 
		* @param new_qpdata data of quadratic program
		*/
		CGQPSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		* 
		*/
		~CGQPSolver();

		void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		void printcontent(ConsoleOutput &output) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printtimer(ConsoleOutput &output) const;
		
		std::string get_name() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(){
	LOG_FUNC_BEGIN

	this->qpdata = NULL;
	
	/* temp vectors */
	g = NULL;
	p = NULL;
	Ap = NULL;
	
	/* iterations counters */
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;	

	/* settings */
	consoleArg.set_option_value("cgqpsolver_maxit", &this->maxit, CGQPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("cgqpsolver_eps", &this->eps, CGQPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("cgqpsolver_debugmode", &this->debugmode, CGQPSOLVER_DEFAULT_DEBUGMODE);
	
	/* timers */
	
	/* function value */
	this->fx = std::numeric_limits<double>::max();

	LOG_FUNC_END
}

template<class VectorBase>
CGQPSolver<VectorBase>::CGQPSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	this->qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();

	/* iterations counters */
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;	

	/* settings */
	consoleArg.set_option_value("cgqpsolver_maxit", &this->maxit, CGQPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("cgqpsolver_eps", &this->eps, CGQPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("cgqpsolver_debugmode", &this->debugmode, CGQPSOLVER_DEFAULT_DEBUGMODE);
	
	/* timers */

	this->fx = std::numeric_limits<double>::max();

	LOG_FUNC_END
}


/* destructor */
template<class VectorBase>
CGQPSolver<VectorBase>::~CGQPSolver(){
	LOG_FUNC_BEGIN

	/* free temp vectors */
	free_temp_vectors();

	LOG_FUNC_END
}

/* prepare temp_vectors */
template<class VectorBase>
void CGQPSolver<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<VectorBase> *pattern = this->qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<VectorBase>(*pattern);
	p = new GeneralVector<VectorBase>(*pattern);
	Ap = new GeneralVector<VectorBase>(*pattern);	
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void CGQPSolver<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(g);
	free(p);
	free(Ap);
	
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	/* print settings */
	if(this->qpdata){
		coutMaster.push();
		this->qpdata->print(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
		
	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void CGQPSolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(this->qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		this->qpdata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void CGQPSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it:          " << this->it_last << std::endl;
	output <<  " - hess mult:   " << this->hessmult_last << std::endl;
	output <<  " - fx:          " << this->fx << std::endl;	
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void CGQPSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - fx:           " << std::setw(25) << this->fx << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void CGQPSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

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

	LOG_FUNC_END
}


template<class VectorBase>
std::string CGQPSolver<VectorBase>::get_name() const {
	std::string return_value = "CGQPSolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}


/* solve the problem */
template<class VectorBase>
void CGQPSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(this->qpdata->get_A());
	pVector b = *(this->qpdata->get_b());
	pVector x0 = *(this->qpdata->get_x0());

	/* pointer to solution */
	pVector x = *(this->qpdata->get_x());

	/* auxiliary vectors */
	pVector g = *(this->g); /* gradient */
	pVector p = *(this->p); /* A-conjugate vector */
	pVector Ap = *(this->Ap); /* A*p */

	x = x0; /* set approximation as initial */

	int it = 0; /* iteration counter */
	int hessmult = 0; /* number of hessian multiplications */
	double normg, alpha, beta, pAp, gg, gg_old;
	
	g = A*x; hessmult += 1; /* compute gradient */
	g -= b;

	p = g; /* initial conjugate direction */

	gg = dot(g,g);

	normg = std::sqrt(gg);

	while(normg > this->eps && it < this->maxit){
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
		
		if(this->debugmode >= 10){
			coutMaster << "it " << it << ": ||g|| = " << normg << std::endl;
		}

		if(this->debugmode >= 100){
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
	if(this->debugmode >= 10){
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

	/* write info to log file */
	LOG_IT(it)
	LOG_FX(this->fx)
		
	LOG_FUNC_END
}

template<class VectorBase>
double CGQPSolver<VectorBase>::get_fx() const {
	if(this->debugmode >= 11) coutMaster << "(CGQPSolver)FUNCTION: get_fx()" << std::endl;
	
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

#ifdef USE_PETSC

/* prepare temp_vectors */
template<>
void CGQPSolver<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	/* I will allocate temp vectors subject to linear term */
	Vec g_vec;
	Vec p_vec;
	Vec Ap_vec;

	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&g_vec) );
	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&p_vec) );
	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&Ap_vec) );

	g = new GeneralVector<PetscVector>(g_vec);
	p = new GeneralVector<PetscVector>(p_vec);
	Ap = new GeneralVector<PetscVector>(Ap_vec);	
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<>
void CGQPSolver<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	Vec g_vec = g->get_vector();
	Vec p_vec = p->get_vector();
	Vec Ap_vec = Ap->get_vector();
	
	TRYCXX( VecDestroy(&g_vec) );
	TRYCXX( VecDestroy(&p_vec) );
	TRYCXX( VecDestroy(&Ap_vec) );

	free(g);
	free(p);
	free(Ap);
	
	LOG_FUNC_END
}

#endif


}
} /* end namespace */

#endif
