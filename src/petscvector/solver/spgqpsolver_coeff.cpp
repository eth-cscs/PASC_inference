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

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();

	/* without cuda */
	TRYCXX( PetscMalloc1(3,&Mdots_val) );
	TRYCXX( PetscMalloc1(3,&(externalcontent->Mdots_vec)) );

	externalcontent->Mdots_vec[0] = d->get_vector();
	externalcontent->Mdots_vec[1] = Ad->get_vector();
	externalcontent->Mdots_vec[2] = g->get_vector();

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
	TRYCXX( PetscFree(externalcontent->Mdots_vec) );
	
	LOG_FUNC_END
}

template<>
void SPGQPSolverC<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const {
	LOG_FUNC_BEGIN

	TRYCXX(VecMDot( externalcontent->Mdots_vec[0], 3, externalcontent->Mdots_vec, Mdots_val) );

	*dd = Mdots_val[0];
	*dAd = Mdots_val[1];
	*gd = Mdots_val[2];

	LOG_FUNC_END
}

template<> 
SPGQPSolverC<PetscVector>::ExternalContent * SPGQPSolverC<PetscVector>::get_externalcontent() const {
	return this->externalcontent;	
}


}
} /* end namespace */

