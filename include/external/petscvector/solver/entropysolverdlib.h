#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H

#include "general/solver/entropysolverdlib.h"

//#include "external/petscvector/solver/generalsolver.h"
//#include "external/petscvector/data/entropydata.h"

#ifdef USE_DLIB

namespace pascinference {
namespace solver {

template<> EntropySolverDlib<PetscVector>::EntropySolverDlib(EntropyData<PetscVector> &new_entropydata);
template<> void EntropySolverDlib<PetscVector>::solve();
template<> void EntropySolverDlib<PetscVector>::compute_moments();
template<> void EntropySolverDlib<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const;

}
} /* end namespace */

#endif
#endif
