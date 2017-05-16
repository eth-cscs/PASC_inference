/** @file qpsolver.h
 *  @brief General Quadratic Programming solver
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_QPSOLVER_H
#define	PASC_QPSOLVER_H

#include "general/solver/generalsolver.h"
#include "general/data/qpdata.h"

#define QPSOLVER_DEFAULT_MAXIT 1000
#define QPSOLVER_DEFAULT_EPS 1e-9
#define QPSOLVER_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class QPSolver
 *  \brief General Quadratic Programming solver
 *
*/
template<class VectorBase>
class QPSolver: public GeneralSolver {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		QPData<VectorBase> *qpdata; /**< data on which the solver operates */
		
		double fx; /**< function value in actual iteration */
		int it_sum; /**< number of all iterations */
		int it_last; /**< number of iterations during last solution */
		int hessmult_sum; /**< number of all Hessian multiplication */
		int hessmult_last; /**< number of Hessian multiplication */

	public:
		/** @brief default constructor
		 * 
		 * Set inner data to null.
		 * 
		 */ 
		QPSolver();
		
		/** @brief constructor from data
		 * 
		 * Set inner pointer to data on which the solver operates to given data.
		 * 
		 * @param new_qpdata input data for solver
		 */ 		
		QPSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		 * 
		 */ 
		~QPSolver();

		virtual void solve();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printstatus(ConsoleOutput &output) const;
		virtual void printstatus(std::ostringstream &output) const;
		virtual void printtimer(ConsoleOutput &output) const;
		virtual void printshort(std::ostringstream &header, std::ostringstream &values) const;
		virtual void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		virtual std::string get_name() const;

		virtual QPData<VectorBase> *get_data() const;
		
		/** @brief return function value in last computed iteration
		* 
		*/		
		virtual double get_fx() const;

		/** @brief return number of iterations during last solve()
		* 
		*/		
		virtual int get_it() const;

		/** @brief return number of Hessian multiplication during last solve()
		* 
		*/		
		virtual int get_hessmult() const;

		ExternalContent *get_externalcontent() const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
QPSolver<VectorBase>::QPSolver(){
	LOG_FUNC_BEGIN
	
	qpdata = NULL;
	
	fx = std::numeric_limits<double>::max();

	consoleArg.set_option_value("qpsolver_maxit", &this->maxit, QPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("qpsolver_eps", &this->eps, QPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("qpsolver_debugmode", &this->debugmode, QPSOLVER_DEFAULT_DEBUGMODE);

	LOG_FUNC_END
}

template<class VectorBase>
QPSolver<VectorBase>::QPSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;
	fx = std::numeric_limits<double>::max();

	consoleArg.set_option_value("qpsolver_maxit", &this->maxit, QPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("qpsolver_eps", &this->eps, QPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("qpsolver_debugmode", &this->debugmode, QPSOLVER_DEFAULT_DEBUGMODE);
			
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
QPSolver<VectorBase>::~QPSolver(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit:       " << this->maxit << std::endl;
	output_global <<  " - eps:         " << this->eps << std::endl;
	output_global <<  " - debugmode:  " << this->debugmode << std::endl;

	/* print data */
	if(qpdata){
		output_global.push();
		qpdata->print(output_global, output_local);
		output_global.pop();
	}
	
	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void QPSolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

/* print status */
template<class VectorBase>
void QPSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* print timer */
template<class VectorBase>
void QPSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string QPSolver<VectorBase>::get_name() const {
	std::string return_value = "QPSolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

/* solve the problem */
template<class VectorBase>
void QPSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN
		
	LOG_FUNC_END
}

template<class VectorBase>
double QPSolver<VectorBase>::get_fx() const {
	return this->fx;
}

template<class VectorBase>
int QPSolver<VectorBase>::get_it() const {
	return this->it_sum;
}

template<class VectorBase>
int QPSolver<VectorBase>::get_hessmult() const {
	return this->hessmult_sum;
}

template<class VectorBase>
QPData<VectorBase> *QPSolver<VectorBase>::get_data() const {
	return this->qpdata;
}


}
} /* end namespace */

#endif
