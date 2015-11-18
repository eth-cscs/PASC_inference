#ifndef QPPROBLEMPROJECTIONSTEP_H
#define	QPPROBLEMPROJECTIONSTEP_H

#include "qpsolver.h"

class QPSolverProjectionstep: public QPSolver {
	protected:
		PetscInt K, N, N_local; /* dimensions of the problem */
		
		Mat Asub; /* block of Hessian matrix */
		Vec *bs; /* blocks of vectors of linear term */
		
		Vec *gs; /* blocks of gradient in x */
		Vec temp,temp2; /* temp vector for computation, same size as gamma_vec[0] */

		PetscScalar stepsize; /* step-size of projection step */

		PetscErrorCode assemble_Asub();
		PetscErrorCode set_bs(Vec *bs);
		PetscErrorCode get_bs(Vec **bs);
		PetscErrorCode compute_gradient();
		PetscErrorCode project();
				
	public:
		QPSolverProjectionstep(Data*, Gamma*, Theta*, PetscScalar);
		PetscErrorCode init();
		PetscErrorCode finalize();
		PetscErrorCode solve();
		PetscErrorCode get_function_value(PetscScalar*);
		PetscErrorCode print(PetscViewer);
		PetscErrorCode correct(PetscScalar);
		
};



#endif		
