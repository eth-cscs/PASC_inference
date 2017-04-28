#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERSPG_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERSPG_H

#include "general/solver/entropysolverspg.h"

//#include "external/petscvector/solver/generalsolver.h"
//#include "external/petscvector/data/entropydata.h"
//#include "external/petscvector/solver/spg_fs.h"
//#include "external/petscvector/algebra/integration/entropyintegration.h"


#ifdef USE_DLIB

namespace pascinference {
namespace solver {

template<> class EntropySolverSPG<PetscVector>::ExternalContent {
	int mynumber;
};

template<> void EntropySolverSPG<PetscVector>::allocate_temp_vectors();
template<> void EntropySolverSPG<PetscVector>::free_temp_vectors();

template<> void EntropySolverSPG<PetscVector>::solve();
template<> void EntropySolverSPG<PetscVector>::compute_moments_data();

template<> void EntropySolverSPG<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;



template<> void EntropySolverSPG<PetscVector>::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec);
template<> double EntropySolverSPG<PetscVector>::compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec);
template<> void EntropySolverSPG<PetscVector>::compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all);


}
} /* end namespace */

#endif
#endif


