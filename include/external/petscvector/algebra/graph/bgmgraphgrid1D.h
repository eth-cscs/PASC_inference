#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPHGRID1D_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/graph/bgmgraphgrid1D.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BGMGraphGrid1D<PetscVector>::ExternalContent : public BGMGraph<PetscVector>::ExternalContent {
	public:
		#ifdef USE_CUDA
			void process_grid_cuda();
		#endif
};

template<> BGMGraphGrid1D<PetscVector>::BGMGraphGrid1D(int width);
template<> void BGMGraphGrid1D<PetscVector>::process_grid();

template<> BGMGraphGrid1D<PetscVector>::ExternalContent * BGMGraphGrid1D<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
