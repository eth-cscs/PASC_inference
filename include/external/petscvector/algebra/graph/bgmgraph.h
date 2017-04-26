#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPH_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPH_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

template<> BGMGraph<PetscVector>::BGMGraph(const double *coordinates_array, int n, int dim);
template<> void BGMGraph<PetscVector>::process(double threshold);
template<> void BGMGraph<PetscVector>::saveVTK(std::string filename) const;

}
} /* end of namespace */

#endif
