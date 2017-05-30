/** @file permonsolver.h
 *  @brief solve QP with solvers implemented in Permon library
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_PERMONSOLVER_H
#define	PASC_PERMONSOLVER_H

#include <iostream>

#include "general/solver/qpsolver.h"
#include "general/data/qpdata.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"

#define PERMONSOLVER_DEFAULT_MAXIT 1000
#define PERMONSOLVER_DEFAULT_EPS 1e-9
#define PERMONSOLVER_USE_UPPERBOUND false
#define PERMONSOLVER_USE_LAMBDAMAX false
#define PERMONSOLVER_DUMP false

namespace pascinference {
namespace solver {

/** \class PermonSolver
 *  \brief Interface with QP solvers implemented in Permon library
 *
 *  For solving QP on closed convex set described by separable simplexes.
*/
template<class VectorBase>
class PermonSolver: public QPSolver<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();


		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		QPData<VectorBase> *qpdata; /**< data on which the solver operates */
		double gP; 					/**< norm of projected gradient */

		/* temporary vectors used during the solution process */
		GeneralVector<VectorBase> *Ad; 		/**< A*d */

		Timer timer_solve; 			/**< total solution time of used algorithm */

		bool dump_or_not;
		bool use_upperbound;		/**< use additional upper bound x<=1 */
		bool use_lambdamax;			/**< provide lambdamax to permon */

	public:
		/** @brief general constructor
		* 
		*/
		PermonSolver();

		/** @brief constructor based on provided data of problem
		* 
		* @param new_qpdata data of quadratic program
		*/
		PermonSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		* 
		*/
		~PermonSolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printtimer(ConsoleOutput &output) const;
		void printshort(std::ostringstream &header, std::ostringstream &values) const;
		void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		std::string get_name() const;

		ExternalContent *get_externalcontent() const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

template<class VectorBase>
void PermonSolver<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("permonsolver_maxit", &this->maxit, PERMONSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("permonsolver_eps", &this->eps, PERMONSOLVER_DEFAULT_EPS);
	
	consoleArg.set_option_value("permonsolver_use_upperbound", &this->use_upperbound, PERMONSOLVER_USE_UPPERBOUND);	
	consoleArg.set_option_value("permonsolver_use_lambdamax", &this->use_lambdamax, PERMONSOLVER_USE_LAMBDAMAX);	
	consoleArg.set_option_value("permonsolver_dump", &this->dump_or_not, PERMONSOLVER_DUMP);	
}


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
PermonSolver<VectorBase>::PermonSolver(){
	LOG_FUNC_BEGIN

	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;
	
	this->fx = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	

	LOG_FUNC_END
}

template<class VectorBase>
PermonSolver<VectorBase>::PermonSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	this->qpdata = &new_qpdata;
	
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
	this->timer_solve.restart();	

	//TODO: in general? Permon is "only" for petsc, this class doesn't make any sence

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
PermonSolver<VectorBase>::~PermonSolver(){
	LOG_FUNC_BEGIN

	/* free temp vectors */
	free_temp_vectors();
	
	LOG_FUNC_END
}

/* prepare temp_vectors */
template<class VectorBase>
void PermonSolver<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<VectorBase> *pattern = this->qpdata->get_b(); /* I will allocate temp vectors subject to linear term */
	Ad = new GeneralVector<VectorBase>(*pattern);	
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void PermonSolver<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(Ad);
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void PermonSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
/*	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	output <<  " - m:          " << m << std::endl;
	output <<  " - gamma:      " << gamma << std::endl;
	output <<  " - sigma1:     " << sigma1 << std::endl;
	output <<  " - sigma2:     " << sigma2 << std::endl;
	output <<  " - alphainit:  " << alphainit << std::endl;
*/	
	/* print data */
	if(this->qpdata){
		coutMaster.push();
		this->qpdata->print(output);
		coutMaster.pop();
	}

	output.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_local <<  " - maxit:      " << this->maxit << std::endl;
	output_local <<  " - eps:        " << this->eps << std::endl;
	output_local <<  " - debugmode: " << this->debugmode << std::endl;

	output_local.synchronize();

	/* print data */
	if(this->qpdata){
		coutMaster.push();
		this->qpdata->print(output_global, output_local);
		coutMaster.pop();
	}

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  " - it: " << std::setw(6) << this->it_last << ", ";
	output <<  "hess mult: " << std::setw(6) << this->hessmult_last << ", ";
	output <<  "fx: " << std::setw(10) << this->fx << ", ";	
//	output <<  "norm(gP): " << std::setw(10) << this->gP << ", ";
	output <<  "used memory: " << std::setw(6) << MemoryCheck::get_virtual() << "%" << std::endl;

	output << " - ";
//	output <<  "t_solve = " << std::setw(10) << this->timer_solve.get_value_last() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(this->qpdata->get_A());
	pVector b = *(this->qpdata->get_b());

	/* pointer to solution */
	pVector x = *(this->qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - fx:           " << std::setw(25) << this->fx << std::endl;
	output <<  "      - fx_control:   " << std::setw(25) << fx_quadratic+fx_linear << ", log: " << std::setw(25) << log(fx_quadratic+fx_linear)/log(10) << std::endl;
	output <<  "      - fx_linear:    " << std::setw(25) << fx_linear << ", log: " << std::setw(25) << log(fx_linear)/log(10) << std::endl;
	output <<  "      - fx_quadratic: " << std::setw(25) << fx_quadratic << ", log: " << std::setw(25) << log(fx_quadratic)/log(10) << std::endl;
//	output <<  "      - norm(gP):     " << std::setw(25) << this->gP << ", log: " << std::setw(25) << log(this->gP)/log(10) << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it all =        " << this->it_sum << std::endl;
	output <<  " - hessmult all =  " << this->hessmult_sum << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =      " << this->timer_solve.get_value_sum() << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(this->qpdata->get_A());
	pVector b = *(this->qpdata->get_b());

	/* pointer to solution */
	pVector x = *(this->qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);
	std::streamsize ss = std::cout.precision();

	values << std::setprecision(17);

	header << "PERMON it, ";
	values << this->it_last << ", ";

	header << "PERMON hessmult, ";
	values << this->hessmult_last << ", ";

	header << "PERMON t all, ";
	values << this->timer_solve.get_value_last() << ", ";

	header << "PERMON fx, ";
	values << this->fx << ", ";

	header << "PERMON fx_linear, ";
	values << fx_linear << ", ";

	header << "PERMON fx_quadratic, ";
	values << fx_quadratic << ", ";

	values << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	header << "PERMON_sum it, ";
	values << this->it_sum << ", ";

	header << "PERMON_sum hessmult, ";
	values << this->hessmult_sum << ", ";

	header << "PERMON_sum t all, ";
	values << this->timer_solve.get_value_sum() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
std::string PermonSolver<VectorBase>::get_name() const {
	std::string return_value = "PermonSolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

/* solve the problem */
template<class VectorBase>
void PermonSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	//TODO: without petsc and permon, this function does not make any sence
	
	LOG_FUNC_END
}



}
} /* end namespace */

#endif
