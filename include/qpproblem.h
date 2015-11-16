#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"

/* include Permon headers */
#include <flloppc.h>
#include <fllopqp.h>
#include <fllopqps.h>
#include <fllopmat.h>

class QPproblem {
		PetscInt K, N, N_local; /* dimensions of the problem */
		
		Mat A; /* Hessian matrix */
		Vec b; /* vetor of linear term */
		
		Mat BE; /* matrix of equality constraints */
		Vec cE; /* vector of equality constraints */
		
		Vec lb; /* lower bound */
		
		Vec x; /* solution vector */
		Vec g; /* gradient in x */
		Vec temp; /* temp vector for computation, same size as x */
		Vec temp2; /* temp vector for computation, same size as cE */

		PetscScalar eps_sqr;
		Mat PBE; /* projector onto Ker B */

	public:

		PetscErrorCode init(PetscInt N, PetscInt N_local, PetscInt K, PetscScalar eps_sqr); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */

		PetscErrorCode assemble_A();
		PetscErrorCode assemble_BE();
		PetscErrorCode assemble_PBE();
		PetscErrorCode assemble_cE();
		PetscErrorCode assemble_lb();

		PetscErrorCode set_b(Vec b);
		PetscErrorCode get_b(Vec *b);
		PetscErrorCode set_x(Vec x);
		PetscErrorCode get_x(Vec *x);

		PetscErrorCode solve_permon();
		PetscErrorCode solve_projection_step();

		PetscErrorCode print(PetscViewer v);

		PetscErrorCode compute_gradient();
		PetscErrorCode get_function_value(PetscScalar *fx);
		PetscErrorCode project();
	
};



#endif
