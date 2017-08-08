#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERNEWTON_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERNEWTON_H

#include "general/solver/entropysolvernewton.h"

/* include integration algorithms */
#include "external/petscvector/algebra/integration/entropyintegration.h"
#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"
//#include "external/petscvector/algebra/integration/entropyintegrationcuba.h"

//#include "external/petscvector/solver/generalsolver.h"
#include "external/petscvector/data/entropydata.h"

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

		Mat H_petsc;						/**< Hessian matrix */
		KSP ksp;							/**< linear solver context */
		PC pc;           					/**< preconditioner context **/

};

template<> EntropySolverNewton<PetscVector>::EntropySolverNewton(EntropyData<PetscVector> &new_entropydata);
template<> EntropySolverNewton<PetscVector>::~EntropySolverNewton();

template<> void EntropySolverNewton<PetscVector>::allocate_temp_vectors();
template<> void EntropySolverNewton<PetscVector>::free_temp_vectors();

template<> void EntropySolverNewton<PetscVector>::solve();

template<> EntropySolverNewton<PetscVector>::ExternalContent * EntropySolverNewton<PetscVector>::get_externalcontent() const;


}
} /* end namespace */

#endif
#endif


