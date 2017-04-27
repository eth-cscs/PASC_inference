#include "external/petscvector/solver/cgqpsolver.h"

namespace pascinference {
namespace solver {

/* prepare temp_vectors */
template<>
void CGQPSolver<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	/* I will allocate temp vectors subject to linear term */
	Vec g_vec;
	Vec p_vec;
	Vec Ap_vec;

	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&g_vec) );
	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&p_vec) );
	TRYCXX( VecDuplicate(this->qpdata->get_b()->get_vector(),&Ap_vec) );

	g = new GeneralVector<PetscVector>(g_vec);
	p = new GeneralVector<PetscVector>(p_vec);
	Ap = new GeneralVector<PetscVector>(Ap_vec);	
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<>
void CGQPSolver<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	Vec g_vec = g->get_vector();
	Vec p_vec = p->get_vector();
	Vec Ap_vec = Ap->get_vector();
	
	TRYCXX( VecDestroy(&g_vec) );
	TRYCXX( VecDestroy(&p_vec) );
	TRYCXX( VecDestroy(&Ap_vec) );

	free(g);
	free(p);
	free(Ap);
	
	LOG_FUNC_END
}


}
} /* end namespace */

