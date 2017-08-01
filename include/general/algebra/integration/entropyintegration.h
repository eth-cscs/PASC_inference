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
#include "general/common/timer.h"

#include "general/data/entropydata.h"

namespace pascinference {
using namespace common;
using namespace data;

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

		EntropyData<VectorBase> *entropydata;
		
		Timer timer;
	public:
		EntropyIntegration(EntropyData<VectorBase> *entropydata, double new_eps);
		~EntropyIntegration();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		double get_eps() const;
		int get_number_of_moments() const;
		int get_xdim() const;
		
		virtual void compute(double *integrals_out, double *lambda, int Km_max = -1);

		double get_time() const;
};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
EntropyIntegration<VectorBase>::EntropyIntegration(EntropyData<VectorBase> *entropydata, double new_eps) {
	LOG_FUNC_BEGIN

	/* set given parameters */
	this->entropydata = entropydata;
	this->eps = new_eps;

	this->timer.restart();
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegration<VectorBase>::~EntropyIntegration(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegration<VectorBase>::compute(double *integrals_out, double *lambda, int Km_max){
	LOG_FUNC_BEGIN
	
	timer.start();
	//TODO: throw error
	timer.stop();
	
	LOG_FUNC_END
}


/* print info about integration */
template<class VectorBase>
void EntropyIntegration<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	output <<  " - number of moments : " << this->entropydata->get_number_of_moments() << std::endl;
	output <<  " - xdim              : " << this->entropydata->get_xdim() << std::endl;
	output <<  " - eps               : " << this->eps << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegration<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	output_global <<  " - number of moments : " << this->entropydata->get_number_of_moments() << std::endl;
	output_global <<  " - xdim              : " << this->entropydata->get_xdim() << std::endl;
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

template<class VectorBase>
int EntropyIntegration<VectorBase>::get_number_of_moments() const {
	return this->entropydata->get_number_of_moments();
}

template<class VectorBase>
int EntropyIntegration<VectorBase>::get_xdim() const {
	return this->entropydata->get_xdim();
}

template<class VectorBase>
double EntropyIntegration<VectorBase>::get_time() const {
	return timer.get_value_sum();
}


}
} /* end namespace */


#endif



