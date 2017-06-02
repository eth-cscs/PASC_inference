/** @file entropyintegration.h
 *  @brief Computes numerical integrals in entropy problem
 *
 *  This integral is computed locally using only "local" MPI process.
 * 
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYINTEGRATION_H
#define	PASC_ENTROPYINTEGRATION_H

#include <string>
#include <iostream>
#include "general/common/consoleoutput.h"
#include "general/common/logging.h"

namespace pascinference {
using namespace common;

namespace algebra {

/** \class EntropyIntegration
 *  \brief Main class providing algorithm for numerical computation of integrals in Entropy problem.
 *
 * \f[
 *  \int\limits_{-1}^{1} x^j e^{-\sum\limits_{k=1}^{m} \lambda_k x^k } ~\textrm{dx}, ~~~ j = 0,\dots,Km
 *	\f] 
 * 
*/
template<class VectorBase>
class EntropyIntegration {
	protected:
		double eps;
	public:
		EntropyIntegration(double new_eps);
		~EntropyIntegration();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		double get_eps() const;
		
		virtual void compute(double *integrals_out, int Km, double *lambda, int Km_max = -1);
};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
EntropyIntegration<VectorBase>::EntropyIntegration(double new_eps) {
	LOG_FUNC_BEGIN

	/* set given parameters */
	this->eps = new_eps;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegration<VectorBase>::~EntropyIntegration(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegration<VectorBase>::compute(double *integrals_out, int Km, double *lambda, int Km_max){
	LOG_FUNC_BEGIN
	
	//TODO
	
	LOG_FUNC_END
}


/* print info about integration */
template<class VectorBase>
void EntropyIntegration<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	output <<  " - eps               : " << this->eps << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegration<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	output_global <<  " - eps               : " << this->eps << std::endl;

	output_global.synchronize();

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyIntegration<VectorBase>::get_name() const {
	return "Entropy-Integration general";
}

template<class VectorBase>
double EntropyIntegration<VectorBase>::get_eps() const {
	return this->eps;
}


}
} /* end namespace */


#endif



