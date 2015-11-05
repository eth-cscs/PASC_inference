#ifndef GAMMA_H
#define	GAMMA_H

class Gamma;

#include "common.h"
#include "data.h"
#include "theta.h"
#include "qpproblem.h"

class Gamma {
		PetscInt dim; /* number of gamma components = K */

		PetscInt global_size; /* global data size with overlap */
		PetscInt local_size; /* local data size with overlap */

		PetscInt local_begin, local_end; /* ownership range */

		PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	
		PetscRandom rnd; /* random numbers generator */

		QPproblem *qpproblem; /* this is qp problem which need to be solved to obtain new gamma */

		/* private functions */
		PetscErrorCode compute_g(Vec g, Data data, Theta theta);
		
	public:
		Vec *gamma_vecs; /* array with data vectors, TODO: should be private */

		PetscErrorCode init(Data, PetscInt, QPproblem*); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */

		PetscErrorCode print(PetscViewer v);

		PetscErrorCode prepare_random();
		PetscErrorCode prepare_uniform();		
		PetscErrorCode prepare_fixed();		

		PetscErrorCode compute(Data data, Theta theta);
		PetscErrorCode projection_step(Data data, Theta theta);
		PetscErrorCode get_objectfunc_value(PetscScalar *value);

		/* GET functions */
		PetscInt get_local_size();
		PetscInt get_local_begin();
		PetscInt get_local_end();
		PetscInt get_global_size();
		PetscInt get_dim();


	
};

#endif

