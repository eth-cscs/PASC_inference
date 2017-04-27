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

template<> void EntropySolverNewton<PetscVector>::allocate_temp_vectors();
template<> void EntropySolverNewton<PetscVector>::free_temp_vectors();

template<> void EntropySolverNewton<PetscVector>::solve();

template<> void EntropySolverNewton<PetscVector>::compute_moments_data();
template<> void EntropySolverNewton<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;


template<> void EntropySolverNewton<PetscVector>::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);
template<> void EntropySolverNewton<PetscVector>::compute_hessian(Vec &integrals_Vec);
template<> double EntropySolverNewton<PetscVector>::compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);
template<> void EntropySolverNewton<PetscVector>::compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all);

}
} /* end namespace */

#endif
#endif


