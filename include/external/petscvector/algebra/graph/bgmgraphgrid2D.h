#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID2D_H

#include "general/algebra/graph/bgmgraphgrid2D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BGMGraphGrid2D<PetscVector>::ExternalContent : public BGMGraph<PetscVector>::ExternalContent {
	public:
		#ifdef USE_CUDA
			void process_grid_cuda(int *neighbor_nmbs, int **neighbor_ids);
		#endif
};

template<> BGMGraphGrid2D<PetscVector>::BGMGraphGrid2D(int width, int height);
template<> void BGMGraphGrid2D<PetscVector>::process_grid();

template<> void BGMGraphGrid2D<PetscVector>::saveVTK_bounding_box(std::string filename, int *bounding_box) const;

template<> BGMGraphGrid2D<PetscVector>::ExternalContent * BGMGraphGrid2D<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
