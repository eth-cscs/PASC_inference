#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"


class QPproblem {
	protected:
		Gamma *gamma;
		Theta *theta;
	public:
		QPproblem(Gamma *, Theta *, PetscScalar);
		virtual PetscErrorCode init();
		virtual PetscErrorCode finalize();
		virtual PetscErrorCode solve();
		virtual PetscErrorCode get_function_value(PetscScalar *fx);
		virtual PetscErrorCode print(PetscViewer v);
		
};



#endif
