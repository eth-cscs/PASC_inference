#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"
#include "theta.h"
#include "data.h"


class QPproblem {
	protected:
		Data *data;
		Gamma *gamma;
		Theta *theta;

		PetscScalar eps_sqr;
	public:
		QPproblem(Data*, Gamma *, Theta *, PetscScalar);
		virtual PetscErrorCode init();
		virtual PetscErrorCode finalize();
		virtual PetscErrorCode solve();
		virtual PetscErrorCode get_function_value(PetscScalar*);
		virtual PetscErrorCode print(PetscViewer);
		virtual PetscErrorCode correct(PetscScalar); /* if L_new > L_old change parameters of solver */
		
};



#endif
