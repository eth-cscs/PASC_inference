#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H

#include "general/solver/entropysolverdlib.h"

//#include "external/petscvector/solver/generalsolver.h"
//#include "external/petscvector/data/entropydata.h"

#ifdef USE_DLIB
#ifdef USE_CUBA

/* include Dlib stuff */
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

/* include Cuba stuff */
#include "cuba.h"

/* Dlib column vector */
typedef dlib::matrix<double,0,1> column_vector;

namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class EntropySolverDlib<PetscVector>::ExternalContent {
	private:
		double integration_eps;
	public:
		class Integrator {
			public:
				int NDIM; //dimensions of integral
				int NCOMP;
				int NVEC;
				double EPSREL;
				double EPSABS;
				int VERBOSE; //log output
				int LAST;
				int SEED;
				int MINEVAL;
				int MAXEVAL;
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

				int comp, nregions, neval, fail;
				cubareal integral[1], error[1], prob[1];

				Integrator();
				~Integrator();

				//four methods of integration implemented in CUBA library,
				//more info at http://www.feynarts.de/cuba/
				double computeVegas();
				double computeSuave();
				double computeDivonne();
				double computeCuhre();

				static int Integrand(const int *ndim, const cubareal xx[],
				const int *ncomp, cubareal ff2[], void *userdata);
		};
		
		class ExtraParameters {
			int k;
			column_vector Mom;
			column_vector LM;
			double L0;
			double eps;
			dlib::matrix<double> D;
			int type; /* type of integrant =0,1,2,3 */
			int order; /* row of D matrix */
			int order2;
    
			ExtraParameters();
			ExtraParameters(int _k, column_vector _Mom, column_vector _LM, double _L0, double _eps, dlib::matrix<double> _D, int _type, int _order);
			~ExtraParameters();
			void Copy(ExtraParameters& _ExtraParameters);			
		};

		ExternalContent(double new_integration_eps);
		double gg(double y, int order, const column_vector& LM);
		double get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps);
		column_vector get_functions_grad(const column_vector& LM, const column_vector& Mom, int k);
		dlib::matrix<double> get_functions_hess(const column_vector& LM, const column_vector& Mom, int k);

		Vec *x_powers_Vecs;
		
};

template<> EntropySolverDlib<PetscVector>::EntropySolverDlib(EntropyData<PetscVector> &new_entropydata);
template<> EntropySolverDlib<PetscVector>::~EntropySolverDlib();

template<> void EntropySolverDlib<PetscVector>::solve();
template<> void EntropySolverDlib<PetscVector>::compute_moments();
template<> void EntropySolverDlib<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;

template<> EntropySolverDlib<PetscVector>::ExternalContent * EntropySolverDlib<PetscVector>::get_externalcontent() const;

}
} /* end namespace */

#endif /* Cuba */
#endif /* Dlib */
#endif
