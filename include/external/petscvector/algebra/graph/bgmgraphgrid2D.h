#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H

#include "general/algebra/graph/bgmgraphgrid2D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

template<> BGMGraphGrid2D<PetscVector>::BGMGraphGrid2D(int width, int height);
template<> void BGMGraphGrid2D<PetscVector>::process_grid();

template<> void BGMGraphGrid2D<PetscVector>::saveVTK_bounding_box(std::string filename, int *bounding_box) const;


}
} /* end of namespace */

#endif
