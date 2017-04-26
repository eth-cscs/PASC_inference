#ifndef PASC_PETSCVECTOR_MULTICGSOLVER_H
#define	PASC_PETSCVECTOR_MULTICGSOLVER_H

#include "general/solver/multicg.h"

//#include "external/petscvector/solver/qpsolver.h"
//#include "external/petscvector/solver/cgqpsolver.h"
//#include "external/petscvector/data/qpdata.h"
//#include "external/petscvector/algebra/matrix/blockdiag.h"
//#include "external/petscvector/algebra/matrix/localdense.h"

namespace pascinference {
namespace solver {

template<> void MultiCGSolver<PetscVector>::solve();

}
} /* end namespace */

#endif
