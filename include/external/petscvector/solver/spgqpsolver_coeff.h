#ifndef PASC_PETSCVECTOR_SPGQPSOLVER_COEFF_H
#define	PASC_PETSCVECTOR_SPGQPSOLVER_COEFF_H

#include "general/solver/spgqpsolver_coeff.h"

//#include "external/petscvector/solver/qpsolver.h"
//#include "external/petscvector/data/qpdata.h"

namespace pascinference {
namespace solver {

template<> void SPGQPSolverC<PetscVector>::allocate_temp_vectors();
template<> void SPGQPSolverC<PetscVector>::free_temp_vectors();

template<> void SPGQPSolverC<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const;


}
} /* end namespace */

#endif
