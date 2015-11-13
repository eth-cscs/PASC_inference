#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"

class QPproblem {
		PetscInt K, N, N_local; /* dimensions of the problem */
		
		Mat A; /* Hessian matrix */
		Vec b; /* vetor of linear term */
		
		Mat BE; /* matrix of equality constraints */
		Vec cE; /* vector of equality constraints */
		
		Vec lb; /* lower bound */
		
		Vec x; /* solution vector */

	public:

		PetscErrorCode init(PetscInt N, PetscInt N_local, PetscInt K); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */

		PetscErrorCode assemble_A();
		PetscErrorCode assemble_BE();
		PetscErrorCode assemble_cE();
		PetscErrorCode assemble_lb();

		PetscErrorCode set_b(Vec b);
		PetscErrorCode get_b(Vec *b);
		PetscErrorCode set_x(Vec x);
		PetscErrorCode get_x(Vec *x);

		PetscErrorCode solve_permon();

		PetscErrorCode print(PetscViewer v);
	
};



#endif
