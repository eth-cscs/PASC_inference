#ifndef PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_PETSC_H
#define	PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_PETSC_H

#include "algebra/feasibleset/simplex_lineqbound.h"

namespace pascinference {
namespace algebra {

template class SimplexFeasibleSet_LinEqBound<PetscVector>;

}
} /* end of namespace */


#endif
