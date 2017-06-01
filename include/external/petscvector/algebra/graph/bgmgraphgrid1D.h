#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H

#include "general/algebra/graph/bgmgraphgrid1D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BGMGraphGrid1D<PetscVector>::ExternalContent : public BGMGraph<PetscVector>::ExternalContent {
	public:
		#ifdef USE_CUDA
			void process_grid_cuda(int *neighbor_nmbs, int **neighbor_ids);
		#endif
};

template<> BGMGraphGrid1D<PetscVector>::BGMGraphGrid1D(int width);
template<> void BGMGraphGrid1D<PetscVector>::process_grid();

template<> BGMGraphGrid1D<PetscVector>::ExternalContent * BGMGraphGrid1D<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
