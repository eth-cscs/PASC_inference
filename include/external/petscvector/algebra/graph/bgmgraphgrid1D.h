#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/graph/bgmgraphgrid1D.h"

namespace pascinference {
namespace algebra {

template<> BGMGraphGrid1D<PetscVector>::BGMGraphGrid1D(int width);

}
} /* end of namespace */

#endif
