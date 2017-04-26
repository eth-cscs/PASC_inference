#ifndef PASC_PETSCVECTOR_SIMPLEX_LINEQBOUNDFEASIBLESET_H
#define	PASC_PETSCVECTOR_SIMPLEX_LINEQBOUNDFEASIBLESET_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"

namespace pascinference {
namespace algebra {

template<> SimplexFeasibleSet_LinEqBound<PetscVector>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K);

template<> Mat SimplexFeasibleSet_LinEqBound<PetscVector>::get_B() const;
template<> Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_c() const;
template<> Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_lb() const;


}
} /* end of namespace */


#endif
