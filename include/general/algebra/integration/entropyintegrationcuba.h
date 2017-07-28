/** @file entropyintegrationcuba.h
 *  @brief Computes numerical integrals using cuba library
 *
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYINTEGRATIONCUBA_H
#define	PASC_ENTROPYINTEGRATIONCUBA_H

#include <string>
#include <iostream>
#include "general/algebra/integration/entropyintegration.h"
#include "general/common/consoleoutput.h"
#include "general/common/logging.h"
#include "general/common/timer.h"

#define ENTROPYINTEGRATIONCUBA_DEFAULT_TYPE 0
#define ENTROPYINTEGRATIONCUBA_DEFAULT_MINEVAL 0
#define ENTROPYINTEGRATIONCUBA_DEFAULT_MAXEVAL 5e4
#define ENTROPYINTEGRATIONCUBA_DEFAULT_NSTART 1e4
#define ENTROPYINTEGRATIONCUBA_DEFAULT_NINCREASE 1e4

#define ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION false

namespace pascinference {
using namespace common;

namespace algebra {

template<class VectorBase>
class EntropyIntegrationCuba : public EntropyIntegration<VectorBase> {
	public:
		class ExternalContent;

		int type;					/**< integration type [0=Vegas,1=Suave,2=Divonne,3=Cuhre] */
		int mineval;				/**< the minimum number of integrand evaluations */
		int maxeval;				/**< the maximum number of integrand evaluations */
		int nstart;					/**< number of integrand evaluations to start with */
		int nincrease;				/**< the increase in number of integrand evaluations */

		bool debug_print_integration;

		void set_settings_from_console();

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

	public:
		EntropyIntegrationCuba(int number_of_moments, int xdim, double new_eps);
		~EntropyIntegrationCuba();

		virtual std::string get_name() const;
		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void compute(double *integrals_out, double *lambda, int Km_max = -1);

		std::string get_integration_type_name(int integration_type) const;		

};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::set_settings_from_console() {

	consoleArg.set_option_value("entropyintegrationcuba_type", &this->type, ENTROPYINTEGRATIONCUBA_DEFAULT_TYPE);
	consoleArg.set_option_value("entropyintegrationcuba_mineval", &this->mineval, ENTROPYINTEGRATIONCUBA_DEFAULT_MINEVAL);
	consoleArg.set_option_value("entropyintegrationcuba_maxeval", &this->maxeval, ENTROPYINTEGRATIONCUBA_DEFAULT_MAXEVAL);
	consoleArg.set_option_value("entropyintegrationcuba_nstart", &this->nstart, ENTROPYINTEGRATIONCUBA_DEFAULT_NSTART);
	consoleArg.set_option_value("entropyintegrationcuba_nincrease", &this->nincrease, ENTROPYINTEGRATIONCUBA_DEFAULT_NINCREASE);
	consoleArg.set_option_value("entropyintegrationcuba_nincrease", &this->nincrease, ENTROPYINTEGRATIONCUBA_DEFAULT_NINCREASE);

	consoleArg.set_option_value("entropyintegrationcuba_debug_print_integration",&debug_print_integration, ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION);
}

/* constructor */
template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::EntropyIntegrationCuba(int number_of_moments, int xdim, double new_eps) : EntropyIntegration<VectorBase>(number_of_moments, xdim, new_eps) {
	LOG_FUNC_BEGIN

	/* load parameters from console */
	set_settings_from_console();

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::~EntropyIntegrationCuba(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::compute(double *integrals_out, double *lambda, int Km_max){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyIntegrationCuba<VectorBase>::get_name() const {
	std::string return_value = "EntropyIntegrationCuba<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
std::string EntropyIntegrationCuba<VectorBase>::get_integration_type_name(int integration_type) const {
	std::string return_string = "undefined";
	switch(integration_type){
		case 0: return_string = "Vegas"; break;
		case 1: return_string = "Suave"; break;
		case 2: return_string = "Divonne"; break;
		case 3: return_string = "Cuhre"; break;
	}
	return return_string;
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	output <<  " - number of moments       : " << this->number_of_moments << std::endl;
	output <<  " - xdim                    : " << this->xdim << std::endl;
	output <<  " - eps                     : " << this->eps << std::endl;

	output <<  " - integration_type        : " << get_integration_type_name(this->type) << std::endl;
	output <<  " - integration_mineval     : " << this->mineval << std::endl;
	output <<  " - integration_maxeval     : " << this->maxeval << std::endl;
	output <<  " - integration_nstart      : " << this->nstart << std::endl;
	output <<  " - integration_nincrease   : " << this->nincrease << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	output_global <<  " - number of moments       : " << this->number_of_moments << std::endl;
	output_global <<  " - xdim                    : " << this->xdim << std::endl;
	output_global <<  " - eps                     : " << this->eps << std::endl;

	output_global <<  " - integration_type        : " << get_integration_type_name(this->type) << std::endl;
	output_global <<  " - integration_mineval     : " << this->mineval << std::endl;
	output_global <<  " - integration_maxeval     : " << this->maxeval << std::endl;
	output_global <<  " - integration_nstart      : " << this->nstart << std::endl;
	output_global <<  " - integration_nincrease   : " << this->nincrease << std::endl;

	output_global.synchronize();

	LOG_FUNC_END
}




}
} /* end namespace */


#endif



