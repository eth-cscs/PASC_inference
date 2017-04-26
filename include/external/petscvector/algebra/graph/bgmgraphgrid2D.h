#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H

#include "general/algebra/graph/bgmgraphgrid2D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

template<> BGMGraphGrid2D<PetscVector>::BGMGraphGrid2D(int width, int height);

}
} /* end of namespace */

#endif
