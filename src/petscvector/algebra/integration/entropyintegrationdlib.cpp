#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"

#ifdef USE_DLIB

namespace pascinference {
namespace algebra {

template<>
void EntropyIntegrationDlib<PetscVector>::compute(double *integrals_arr, double *lambda_arr, int Km_max) {
	LOG_FUNC_BEGIN
	
	if(Km_max < 0){
		Km_max = this->entropydata->get_number_of_moments()-1;
	}
	
	column_vector lambda_Dlib(this->entropydata->get_number_of_moments()-1);

	/* from arr to Dlib-vec */
	for(int km=0;this->entropydata->get_number_of_moments()-1;km++){
		lambda_Dlib(km) = lambda_arr[km];
	}

	/* compute integrals */
    for(int km = 0; km<Km_max;km++){
		auto mom_function = [&](double x)->double { return this->externalcontent->gg(x, km, lambda_Dlib);};
		integrals_arr[km] = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->eps);
	}
	
	LOG_FUNC_END
}

double EntropyIntegrationDlib<PetscVector>::ExternalContent::gg(double y, int order, column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    return pow(y,order)*(exp(-trans(LM)*z));
}


}
}

#endif
