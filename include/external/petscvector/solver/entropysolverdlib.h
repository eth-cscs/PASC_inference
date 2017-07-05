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
typedef dlib::matrix<double,1,0> row_vector;

namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class EntropySolverDlib<PetscVector>::ExternalContent {
	private:
		double integration_eps;
		int integration_type;

		column_vector cLM;
		column_vector cgrad;
		double cF;
		
		dlib::matrix<double> chess;

		bool debug_print_content;
		bool debug_print_integration;

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
				
				Integrator(int integration_type, int ndim, int ncomp, bool debug_print_integration = false);
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
				
				Timer timer; /**< total integration time */
				double get_time() const;
		};
		
		class ExtraParameters {
			public:
				column_vector Mom;
				column_vector LM;
				double eps;
				dlib::matrix<double> D;

				ExtraParameters();
				ExtraParameters(column_vector _Mom, column_vector _LM, dlib::matrix<double> _D, double _eps);
				~ExtraParameters();
				void Copy(ExtraParameters& _ExtraParameters);			
		};

		ExternalContent(double new_integration_eps, int integration_type=0, bool debug_print_content = false, bool debug_print_integration = false);
		double gg(double y, int order, const column_vector& LM);
		double get_functions_obj(const column_vector& _LM, const column_vector& _Mom, double eps, int k, const dlib::matrix<double>& mom_powers);
		column_vector get_functions_grad(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers);
		dlib::matrix<double> get_functions_hess(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers);

		double get_F() const;

		Vec *x_powers_Vecs;
		double *Fs; /**< value of F for all clusters */
		double integration_time; /**< total integration time */
};

template<> EntropySolverDlib<PetscVector>::EntropySolverDlib(EntropyData<PetscVector> &new_entropydata);
template<> EntropySolverDlib<PetscVector>::~EntropySolverDlib();

template<> void EntropySolverDlib<PetscVector>::solve();
template<> void EntropySolverDlib<PetscVector>::compute_moments();
template<> void EntropySolverDlib<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;

template<> double EntropySolverDlib<PetscVector>::get_integration_time() const;

template<> EntropySolverDlib<PetscVector>::ExternalContent * EntropySolverDlib<PetscVector>::get_externalcontent() const;

}
} /* end namespace */

#endif /* Cuba */
#endif /* Dlib */
#endif
