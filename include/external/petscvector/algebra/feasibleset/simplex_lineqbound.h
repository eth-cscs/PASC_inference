#ifndef PASC_PETSCVECTOR_SIMPLEX_LINEQBOUNDFEASIBLESET_H
#define	PASC_PETSCVECTOR_SIMPLEX_LINEQBOUNDFEASIBLESET_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class SimplexFeasibleSet_LinEqBound<PetscVector>::ExternalContent {
	public:
		Mat B; /**< matrix of equality constraints Bx=c */
		Vec c; /**< vector of equality constraints Bx=c */
		Vec lb; /**< vector of lower bounds */
};

template<> SimplexFeasibleSet_LinEqBound<PetscVector>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K);
template<> void SimplexFeasibleSet_LinEqBound<PetscVector>::project(GeneralVector<PetscVector> &x);

template<> SimplexFeasibleSet_LinEqBound<PetscVector>::ExternalContent * SimplexFeasibleSet_LinEqBound<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */


#endif
