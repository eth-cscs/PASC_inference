#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID3D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID3D_H

#include "general/algebra/graph/bgmgraphgrid3D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BGMGraphGrid3D<PetscVector>::ExternalContent : public BGMGraph<PetscVector>::ExternalContent {
	public:
		#ifdef USE_CUDA
			void process_grid_cuda(int *neighbor_nmbs, int **neighbor_ids);
		#endif
};

template<> BGMGraphGrid3D<PetscVector>::BGMGraphGrid3D(int x_size, int y_size, int z_size);
template<> void BGMGraphGrid3D<PetscVector>::process_grid();

template<> void BGMGraphGrid3D<PetscVector>::saveVTK_bounding_box(std::string filename, int *bounding_box) const;

template<> BGMGraphGrid3D<PetscVector>::ExternalContent * BGMGraphGrid3D<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
