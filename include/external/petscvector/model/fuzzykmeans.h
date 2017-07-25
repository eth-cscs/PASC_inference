#ifndef PASC_PETSCVECTOR_FUZZYKMEANSMODEL_H
#define PASC_PETSCVECTOR_FUZZYKMEANSMODEL_H

#include "general/model/fuzzykmeans.h"
#include "external/petscvector/common/common.h"

/* gamma, theta problem */
//#include "external/petscvector/solver/simplesolver.h"
//#include "external/petscvector/data/simpledata.h"
#include "external/petscvector/data/tsdata.h"


namespace pascinference {
namespace model {

template<> FuzzyKmeansModel<PetscVector>::FuzzyKmeansModel(TSData<PetscVector> &new_tsdata, double fuzzifier, Fem<PetscVector> *new_fem);
template<> void FuzzyKmeansModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

template<> void FuzzyKmeansModel<PetscVector>::gammasolver_initialize(GeneralSolver **gammasolver);
template<> void FuzzyKmeansModel<PetscVector>::gammasolver_updatebeforesolve(GeneralSolver *gammasolver);
template<> void FuzzyKmeansModel<PetscVector>::gammasolver_updateaftersolve(GeneralSolver *gammasolver);

template<> void FuzzyKmeansModel<PetscVector>::thetasolver_initialize(GeneralSolver **thetasolver);
template<> void FuzzyKmeansModel<PetscVector>::thetasolver_updatebeforesolve(GeneralSolver *thetasolver);
template<> void FuzzyKmeansModel<PetscVector>::thetasolver_updateaftersolve(GeneralSolver *thetasolver);

template<> double FuzzyKmeansModel<PetscVector>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver);
template<> void FuzzyKmeansModel<PetscVector>::compute_residuum();

}
} /* end namespace */

#endif
