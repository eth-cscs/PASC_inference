#include "external/petscvector/solver/spgqpsolver_coeff.h"

namespace pascinference {
namespace solver {

/* prepare temp_vectors */
template<>
void SPGQPSolverC<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<PetscVector> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<PetscVector>(*pattern);
	d = new GeneralVector<PetscVector>(*pattern);
	Ad = new GeneralVector<PetscVector>(*pattern);	
	temp = new GeneralVector<PetscVector>(*pattern);	

	/* without cuda */
	TRYCXX( PetscMalloc1(3,&Mdots_val) );
	TRYCXX( PetscMalloc1(3,&Mdots_vec) );

	Mdots_vec[0] = d->get_vector();
	Mdots_vec[1] = Ad->get_vector();
	Mdots_vec[2] = g->get_vector();

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<>
void SPGQPSolverC<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(g);
	free(d);
	free(Ad);
	free(temp);
	
	TRYCXX( PetscFree(Mdots_val) );
	TRYCXX( PetscFree(Mdots_vec) );
	
	LOG_FUNC_END
}

template<>
void SPGQPSolverC<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const {
	LOG_FUNC_BEGIN

	TRYCXX(VecMDot( Mdots_vec[0], 3, Mdots_vec, Mdots_val) );

	*dd = Mdots_val[0];
	*dAd = Mdots_val[1];
	*gd = Mdots_val[2];

	LOG_FUNC_END
}


}
} /* end namespace */

