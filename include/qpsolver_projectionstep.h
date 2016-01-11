#ifndef QPPROBLEMPROJECTIONSTEP_H
#define	QPPROBLEMPROJECTIONSTEP_H

#include "qpsolver.h"

class QPSolverProjectionstep: public QPSolver {
	protected:
		int K, N, N_local; /* dimensions of the problem */
		
		Mat Asub; /* block of Hessian matrix */
		Vec *bs; /* blocks of vectors of linear term */
		
		Vec *gs; /* blocks of gradient in x */
		Vec temp,temp2; /* temp vector for computation, same size as gamma_vec[0] */

		PetscScalar stepsize; /* step-size of projection step */

		void assemble_Asub();
		void set_bs(Vec *bs);
		void get_bs(Vec **bs);
		void compute_gradient();
		void project();
				
	public:
		QPSolverProjectionstep(Data*, Gamma*, Theta*, PetscScalar);
		void init();
		void finalize();
		void solve();
		void get_function_value(PetscScalar*);
		void print(PetscViewer);
		void correct(PetscScalar);
		
};



#endif		
