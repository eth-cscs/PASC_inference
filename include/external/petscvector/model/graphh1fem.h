#ifndef PASC_PETSCVECTOR_GRAPHH1FEMMODEL_H
#define PASC_PETSCVECTOR_GRAPHH1FEMMODEL_H

#include "general/model/graphh1fem.h"

#include "external/petscvector/common/common.h"

/* gamma problem */
//#include "external/petscvector/algebra/matrix/blockgraphfree.h" // TODO: implement?
#include "external/petscvector/algebra/matrix/blockgraphsparse.h"

#include "external/petscvector/algebra/feasibleset/simplex_local.h"
#include "external/petscvector/solver/spgqpsolver.h"
//#include "external/petscvector/solver/spgqpsolver_coeff.h"
//#include "external/petscvector/solver/taosolver.h"
//#include "external/petscvector/data/qpdata.h"

//TODO: !!!
#ifdef USE_PERMON
	#include "external/petscvector/solver/permonsolver.h"
	#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"
#endif

/* theta problem */
//#include "external/petscvector/solver/simplesolver.h"
//#include "external/petscvector/data/simpledata.h"
#include "external/petscvector/data/tsdata.h"


namespace pascinference {
namespace model {

template<> GraphH1FEMModel<PetscVector>::GraphH1FEMModel(TSData<PetscVector> &new_tsdata, double epssqr, Fem<PetscVector> *new_fem, bool usethetainpenalty);
template<> void GraphH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

template<> void GraphH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver);
template<> void GraphH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver);

template<> void GraphH1FEMModel<PetscVector>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver);
template<> void GraphH1FEMModel<PetscVector>::updateaftersolve_gammasolver(GeneralSolver *gammasolver);

template<> void GraphH1FEMModel<PetscVector>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver);
template<> void GraphH1FEMModel<PetscVector>::updateaftersolve_thetasolver(GeneralSolver *thetasolver);

}
} /* end namespace */

#endif
