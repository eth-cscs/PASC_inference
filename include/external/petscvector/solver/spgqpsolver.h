#ifndef PASC_PETSCVECTOR_SPGQPSOLVER_H
#define	PASC_PETSCVECTOR_SPGQPSOLVER_H

#include "general/solver/spgqpsolver.h"

#include "external/petscvector/common/common.h"
//#include "external/petscvector/solver/qpsolver.h"
//#include "external/petscvector/data/qpdata.h"
//#include "external/petscvector/solver/spg_fs.h"
#include "external/petscvector/algebra/matrix/blockgraphsparse.h"

#ifdef USE_CUDA
 #include "petsccuda.h"											/* VecCUDAGetArrayReadWrite */
 #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>		/* VecCUDACopyToGPU */
#endif


namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class SPGQPSolver<PetscVector>::ExternalContent {
	public:
		Vec *Mdots_vec; /**< for manipulation with mdot */
};

template<> std::string SPGQPSolver<PetscVector>::get_name() const;
template<> void SPGQPSolver<PetscVector>::allocate_temp_vectors();
template<> void SPGQPSolver<PetscVector>::free_temp_vectors();
template<> void SPGQPSolver<PetscVector>::solve();
template<> double SPGQPSolver<PetscVector>::get_fx() const;
template<> void SPGQPSolver<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const;

template<> SPGQPSolver<PetscVector>::ExternalContent * SPGQPSolver<PetscVector>::get_externalcontent() const;

}
} /* end namespace */


#endif
