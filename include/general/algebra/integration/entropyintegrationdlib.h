/** @file entropyintegrationdlib.h
 *  @brief Computes numerical integrals using dlib library
 *
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYINTEGRATIONDLIB_H
#define	PASC_ENTROPYINTEGRATIONDLIB_H

#include <string>
#include <iostream>
#include "general/algebra/integration/entropyintegration.h"
#include "general/common/consoleoutput.h"
#include "general/common/logging.h"
#include "general/common/timer.h"

namespace pascinference {
using namespace common;

namespace algebra {

template<class VectorBase>
class EntropyIntegrationDlib : public EntropyIntegration<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

	public:
		EntropyIntegrationDlib(int number_of_moments, double new_eps);
		~EntropyIntegrationDlib();

		virtual std::string get_name() const;

		virtual void compute(double *integrals_out, double *lambda, int Km_max = -1);
};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
EntropyIntegrationDlib<VectorBase>::EntropyIntegrationDlib(int number_of_moments, double new_eps) : EntropyIntegration<VectorBase>(number_of_moments, 1, NULL, new_eps) { /* here xdim=1, matrix_D_arr=NULL */
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegrationDlib<VectorBase>::~EntropyIntegrationDlib(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegrationDlib<VectorBase>::compute(double *integrals_out, double *lambda, int Km_max){
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyIntegrationDlib<VectorBase>::get_name() const {
	std::string return_value = "EntropyIntegrationDlib<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}


}
} /* end namespace */


#endif



