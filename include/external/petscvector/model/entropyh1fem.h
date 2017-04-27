#ifndef PASC_PETSCVECTOR_ENTROPYH1FEMMODEL_H
#define PASC_PETSCVECTOR_ENTROPYH1FEMMODEL_H

#include "general/model/entropyh1fem.h"

//#include "external/petscvector/model/tsmodel.h"

/* gamma problem */
#include "external/petscvector/algebra/matrix/blockgraphsparse.h"
//#include "external/petscvector/data/qpdata.h"
#include "external/petscvector/algebra/feasibleset/simplex_local.h"
#include "external/petscvector/solver/spgqpsolver.h"
//#include "external/petscvector/solver/spgqpsolver_coeff.h"
//#include "external/petscvector/solver/taosolver.h"

//TODO: !!!
#ifdef USE_PERMON
	#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"
	#include "external/petscvector/solver/permonsolver.h"
#endif

/* theta problem */
//#include "external/petscvector/data/entropydata.h"
#include "external/petscvector/solver/entropysolverdlib.h"
#include "external/petscvector/solver/entropysolverspg.h"
#include "external/petscvector/solver/entropysolvernewton.h"
#include "external/petscvector/data/tsdata.h"

namespace pascinference {
namespace model {

template<> EntropyH1FEMModel<PetscVector>::EntropyH1FEMModel(TSData<PetscVector> &new_tsdata, int Km, double epssqr);
template<> void EntropyH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

template<> void EntropyH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver);
template<> void EntropyH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver);

template<> void EntropyH1FEMModel<PetscVector>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver);
template<> void EntropyH1FEMModel<PetscVector>::updateaftersolve_gammasolver(GeneralSolver *gammasolver);

template<> void EntropyH1FEMModel<PetscVector>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver);
template<> void EntropyH1FEMModel<PetscVector>::updateaftersolve_thetasolver(GeneralSolver *thetasolver);


}
} /* end namespace */

#endif
