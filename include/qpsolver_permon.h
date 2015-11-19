#ifndef QPPROBLEMPERMON_H
#define	QPPROBLEMPERMON_H

#include "qpsolver.h"

/* include Permon headers */
#include <flloppc.h>
#include <fllopqp.h>
#include <fllopqps.h>
#include <fllopmat.h>


class QPSolverPermon: public QPSolver {
	protected:
		PetscInt K, N, N_local; /* dimensions of the problem */
		
		Mat A; /* Hessian matrix */
		Vec b; /* vector of linear term */
		
		Mat R; /* Im R = Ker A */
		
		Mat BE; /* matrix of equality constraints */
		Vec cE; /* vector of equality constraints */
		
		Vec lb; /* lower bound */
		
		Vec x; /* solution vector */
		Vec g; /* gradient in x */
		Vec temp; /* temp vector for computation, same size as x */
		Vec temp2; /* temp vector for computation, same size as cE */

		PetscScalar rtol;
		
		PetscErrorCode assemble_A();
		PetscErrorCode assemble_R();
		PetscErrorCode assemble_BE();
		PetscErrorCode assemble_cE();
		PetscErrorCode assemble_lb();
		PetscErrorCode set_b(Vec b);
		PetscErrorCode get_b(Vec *b);
		PetscErrorCode set_x(Vec x);
		PetscErrorCode get_x(Vec *x);
		PetscErrorCode compute_gradient();
		PetscErrorCode project();
				
	public:
		QPSolverPermon(Data*, Gamma*, Theta*, PetscScalar);
		PetscErrorCode init();
		PetscErrorCode finalize();
		PetscErrorCode solve();
		PetscErrorCode get_function_value(PetscScalar *fx);
		PetscErrorCode print(PetscViewer v);
		PetscErrorCode correct(PetscScalar);

		
		
};



#endif		
