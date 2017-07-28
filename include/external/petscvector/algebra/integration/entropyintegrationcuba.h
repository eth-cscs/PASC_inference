#ifndef PASC_PETSCVECTOR_ENTROPYINTEGRATIONCUBA_H
#define	PASC_PETSCVECTOR_ENTROPYINTEGRATIONCUBA_H

#include "general/algebra/integration/entropyintegrationcuba.h"
#include "external/petscvector/algebra/vector/generalvector.h"

//#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"


#ifdef USE_CUBA
/* if we are not using CUBA, then this class does not make any sence */

/* include Cuba stuff */
#include "cuba.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class EntropyIntegrationCuba<PetscVector>::ExternalContent {
	public:
		class Integrator {
			public:
				int NDIM; 		/**< number of dimensions of integral */
				int NCOMP;		/**< number of components of the integrand */
				int NVEC;	
				double EPSREL;	/**< requested relative accuracy */
				double EPSABS;	/**< requested absolute accuracy */
				int VERBOSE; //log output
				int LAST;
				int SEED;
				int MINEVAL;	/**< minimum number of integrand evaluations */
				int MAXEVAL;	/**< maximum number of integrand evaluations allowed */
				int NSTART;
				int NINCREASE;
				int NBATCH;
				int GRIDNO;
				char* STATEFILE = NULL;
				void* SPIN = NULL;
				int NNEW;
				int NMIN;
				double FLATNESS;
				void* USERDATA = NULL; //this is to pass extra parameters to integral

				int KEY1;
				int KEY2;
				int KEY3;
				int MAXPASS;
				double BORDER;
				double MAXCHISQ;
				double MINDEVIATION;
				int NGIVEN;
				int LDXGIVEN;
				int NEXTRA;
				int KEY;

				int integration_type;

				int comp, nregions, neval, fail;

				cubareal *integral;
				cubareal *error;
				cubareal *prob;
				
				bool debug_print_integration;
				
				Integrator(int integration_type, int ndim, int ncomp, int integration_mineval, int integration_maxeval, int integration_nstart, int integration_nincrease, bool debug_print_integration = ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION);
				~Integrator();

				//four methods of integration implemented in CUBA library,
				//more info at http://www.feynarts.de/cuba/
				void computeVegas();
				void computeSuave();
				void computeDivonne();
				void computeCuhre();
				cubareal *compute();

				static int Integrand(const int *ndim, const cubareal xx[],
				const int *ncomp, cubareal ff2[], void *userdata);
		};
		
		class ExtraParameters {
			public:
				double *Mom;
				double *LM;
				double eps;
				double *D_matrix;

				ExtraParameters();
				ExtraParameters(double *_Mom, double *_LM, double *_D, double _eps);
				~ExtraParameters();
				void Copy(ExtraParameters& _ExtraParameters);			
		};

};


template<> void EntropyIntegrationCuba<PetscVector>::compute(double *integrals_out, double *lambda, int Km_max);

}
} /* end of namespace */

#endif

#endif
