#ifndef PASC_PETSCVECTOR_GENERALMATRIX_H
#define	PASC_PETSCVECTOR_GENERALMATRIX_H

#include "general/algebra/vector/generalvector.h"
#include "general/algebra/matrix/generalmatrix.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class GeneralMatrix<PetscVector>::ExternalContent {
	public:
		Mat A_petsc; /**< internal PETSc matrix */
};

}
} /* end of namespace */


#endif
