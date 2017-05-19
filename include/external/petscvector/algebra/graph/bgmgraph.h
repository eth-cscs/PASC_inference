#ifndef PASC_PETSCVECTOR_COMMON_BGMGRAPH_H
#define	PASC_PETSCVECTOR_COMMON_BGMGRAPH_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class BGMGraph<PetscVector>::ExternalContent {
	public:
		int n;

		#ifdef USE_CUDA
			int *neighbor_nmbs_gpu; /**< copy of values on GPU */
			int **neighbor_ids_cpugpu; /**< pointers to GPU arrays on CPU */
			int **neighbor_ids_gpu; /**< copy of values on GPU */

			void cuda_destroy();
			void cuda_process();
		#endif
		
};

template<> BGMGraph<PetscVector>::BGMGraph(const double *coordinates_array, int n, int dim);
template<> BGMGraph<PetscVector>::~BGMGraph();

template<> int *BGMGraph<PetscVector>::get_neighbor_nmbs_gpu() const;
template<> int **BGMGraph<PetscVector>::get_neighbor_ids_gpu() const;

template<> void BGMGraph<PetscVector>::process(double threshold);
template<> void BGMGraph<PetscVector>::saveVTK(std::string filename) const;

template<> BGMGraph<PetscVector>::ExternalContent * BGMGraph<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
