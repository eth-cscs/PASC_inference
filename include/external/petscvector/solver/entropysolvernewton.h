#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERNEWTON_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERNEWTON_H

#include "general/solver/entropysolvernewton.h"

//#include "external/petscvector/solver/generalsolver.h"
//#include "external/petscvector/data/entropydata.h"
//#include "external/petscvector/algebra/integration/entropyintegration.h"

/* petsc stuff */
#include <petscksp.h>

#ifdef USE_DLIB

namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class EntropySolverNewton<PetscVector>::ExternalContent {
	public:
		void compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);					/**< g = nabla_lambda f(lambda, moments) */
		void compute_hessian(Vec &integrals_Vec);													/**< Hessian matrix */
		double compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);		/**< compute function value from already computed integrals and moments */
		void compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, EntropyIntegration<PetscVector> *entropyintegration, bool compute_all);				/**< compute integrals int x^{0,..,Km} exp(-dot(lambda,x^{1,..,Km})) */

		Mat H_petsc;						/**< Hessian matrix */
		KSP ksp;							/**< linear solver context */
		PC pc;           					/**< preconditioner context **/

		/* functions for Dlib */
		#ifdef USE_DLIB
			double gg(double y, int order, column_vector& LM);
		#endif

};

template<> EntropySolverNewton<PetscVector>::EntropySolverNewton(EntropyData<PetscVector> &new_entropydata);

template<> void EntropySolverNewton<PetscVector>::allocate_temp_vectors();
template<> void EntropySolverNewton<PetscVector>::free_temp_vectors();

template<> void EntropySolverNewton<PetscVector>::solve();

template<> void EntropySolverNewton<PetscVector>::compute_moments_data();
template<> void EntropySolverNewton<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;

template<> EntropySolverNewton<PetscVector>::ExternalContent * EntropySolverNewton<PetscVector>::get_externalcontent() const;


}
} /* end namespace */

#endif
#endif


