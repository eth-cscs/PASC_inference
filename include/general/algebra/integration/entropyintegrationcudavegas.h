/** @file entropyintegrationcudavegas.h
 *  @brief Computes numerical integrals using cuda implementation of Vegas algorithm
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_ENTROPYINTEGRATIONCUDAVEGAS_H
#define	PASC_ENTROPYINTEGRATIONCUDAVEGAS_H

#include <string>
#include <iostream>
#include "general/algebra/integration/entropyintegration.h"
#include "general/common/consoleoutput.h"
#include "general/common/logging.h"
#include "general/common/timer.h"

#define ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_NCALL 500000
#define ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_ITMX 10
#define ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_NBLOCKSIZE 256
#define ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_DEBUG_PRINT_INTEGRATION false
#define ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_DEBUG_PRINT_INTEGRATION_INNER false

namespace pascinference {
using namespace common;

namespace algebra {

template<class VectorBase>
class EntropyIntegrationCudaVegas : public EntropyIntegration<VectorBase> {
	public:
		class ExternalContent;

		int ncall;
		int itmx;
		int nBlockSize;
		bool debug_print_integration;
		bool debug_print_integration_inner;


	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		void set_settings_from_console();

	public:
		EntropyIntegrationCudaVegas(EntropyData<VectorBase> *entropydata, double new_eps);
		~EntropyIntegrationCudaVegas();

		virtual std::string get_name() const;
		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void compute(double *integrals_out, double *lambda, int Km_max = -1);

};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

template<class VectorBase>
void EntropyIntegrationCudaVegas<VectorBase>::set_settings_from_console() {

	consoleArg.set_option_value("entropyintegrationcudavegas_ncall",&ncall, ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_NCALL);
	consoleArg.set_option_value("entropyintegrationcudavegas_itmx",&itmx, ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_ITMX);
	consoleArg.set_option_value("entropyintegrationcudavegas_nblocksize",&nBlockSize, ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_NBLOCKSIZE);
	consoleArg.set_option_value("entropyintegrationcudavegas_debug_print_integration",&debug_print_integration, ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_DEBUG_PRINT_INTEGRATION);
	consoleArg.set_option_value("entropyintegrationcudavegas_debug_print_integration_inner",&debug_print_integration_inner, ENTROPYINTEGRATIONCUDAVEGAS_DEFAULT_DEBUG_PRINT_INTEGRATION_INNER);

}

/* constructor */
template<class VectorBase>
EntropyIntegrationCudaVegas<VectorBase>::EntropyIntegrationCudaVegas(EntropyData<VectorBase> *entropydata, double new_eps) : EntropyIntegration<VectorBase>(entropydata, new_eps) {
	LOG_FUNC_BEGIN

	/* load parameters from console */
	set_settings_from_console();

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegrationCudaVegas<VectorBase>::~EntropyIntegrationCudaVegas(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyIntegrationCudaVegas<VectorBase>::get_name() const {
	std::string return_value = "EntropyIntegrationCudaVegas<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCudaVegas<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;

	output <<  " - number of moments        : " << this->entropydata->get_number_of_moments() << std::endl;
	output <<  " - xdim                     : " << this->entropydata->get_xdim() << std::endl;
	output <<  " - ncall                    : " << this->ncall << std::endl;
	output <<  " - itmx                     : " << this->itmx << std::endl;
	output <<  " - nblocksize               : " << this->nBlockSize << std::endl;
	output <<  " - eps                      : " << this->eps << std::endl;

	output.synchronize();

	LOG_FUNC_END
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCudaVegas<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	output_global <<  " - number of moments       : " << this->entropydata->get_number_of_moments() << std::endl;
	output_global <<  " - xdim                    : " << this->entropydata->get_xdim() << std::endl;
	output_global <<  " - ncall                   : " << this->ncall << std::endl;
	output_global <<  " - itmx                    : " << this->itmx << std::endl;
	output_global <<  " - nblocksize              : " << this->nBlockSize << std::endl;
	output_global <<  " - eps                     : " << this->eps << std::endl;

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegrationCudaVegas<VectorBase>::compute(double *integrals_out, double *lambda, int Km_max) {
	LOG_FUNC_BEGIN

	this->timer.start();

	//TODO

/*    if(debug_print_integration){
        coutAll << "lambda:    " << print_array(lambda, this->entropydata->get_number_of_moments()-1) << std::endl;
        coutAll << "integrals: " << print_array(computed_integrals, Km_max) << std::endl;
        coutAll.synchronize();
    }
*/
	this->timer.stop();

	LOG_FUNC_END
}


}
} /* end namespace */


#endif



