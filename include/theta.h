#ifndef THETA_H
#define	THETA_H

class Theta;

#include "common.h"
#include "data.h"
#include "gamma.h"

class Theta {
		PetscInt dim_data; /* number of data components = n */
		PetscInt dim_gamma; /* number of gamma components = K */

	public:
		PetscScalar *theta_arr; /* array with theta, should be same in every proc */

		PetscErrorCode init(Data data, Gamma gamma); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */
		PetscErrorCode print(PetscViewer v);
		PetscErrorCode compute(Data data, Gamma gamma);
		
		
		/* GET functions */
		PetscInt get_dim_data();
		PetscInt get_dim_gamma();

	
};

#endif
