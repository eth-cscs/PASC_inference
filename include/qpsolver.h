#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"
#include "theta.h"
#include "data.h"


class QPSolver {
	protected:
		Data *data;
		Gamma *gamma;
		Theta *theta;

		PetscScalar eps_sqr;
	public:
		QPSolver(Data*, Gamma *, Theta *, PetscScalar);
		virtual void init();
		virtual void finalize();
		virtual void solve();
		virtual void get_function_value(PetscScalar*);
		virtual void print(PetscViewer);
		virtual void correct(PetscScalar); /* if L_new > L_old change parameters of solver */
		
};



#endif
