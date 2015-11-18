#ifndef GAMMA_H
#define	GAMMA_H

class Gamma;
class QPSolver;

#include "common.h"
#include "data.h"
#include "theta.h"
#include "qpsolver.h"

class Gamma {
		PetscInt dim; /* number of gamma components = K */

		PetscInt global_size; /* global data size with overlap */
		PetscInt local_size; /* local data size with overlap */

		PetscInt local_begin, local_end; /* ownership range */

		PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	
		PetscRandom rnd; /* random numbers generator */
		
	public:
		Vec *gamma_vecs; /* array with data vectors, TODO: should be private */

		PetscErrorCode init(Data, PetscInt); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */

		PetscErrorCode print(PetscViewer v);

		PetscErrorCode prepare_random();
		PetscErrorCode prepare_uniform();		
		PetscErrorCode prepare_fixed();		

		PetscErrorCode compute(QPSolver *qp_solver, Data data, Theta theta);

		/* GET functions */
		PetscInt get_local_size();
		PetscInt get_local_begin();
		PetscInt get_local_end();
		PetscInt get_global_size();
		PetscInt get_dim();

		// TODO: this should be somewhere else, Model?
		PetscErrorCode compute_g(Vec g, Data *data, Theta *theta);
		PetscErrorCode compute_gk(Vec g, Data *data, Theta *theta, PetscInt k);
	
};

#endif

