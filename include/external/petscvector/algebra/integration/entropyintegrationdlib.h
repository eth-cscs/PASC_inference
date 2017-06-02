#ifndef PASC_PETSCVECTOR_ENTROPYINTEGRATIONDLIB_H
#define	PASC_PETSCVECTOR_ENTROPYINTEGRATIONDLIB_H

#include "general/algebra/integration/entropyintegrationdlib.h"
#include "external/petscvector/algebra/vector/generalvector.h"

//#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"


#ifdef USE_DLIB
/* if we are not using DLIB, then this class does not make any sence */

#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

/* Dlib column vector */
typedef dlib::matrix<double,0,1> column_vector;

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class EntropyIntegrationDlib<PetscVector>::ExternalContent {
	public:
		double gg(double y, int order, column_vector& LM);

};


template<> void EntropyIntegrationDlib<PetscVector>::compute(double *integrals_out, int Km, double *lambda, int Km_max);

}
} /* end of namespace */

#endif

#endif
