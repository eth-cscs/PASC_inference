#ifndef PASCSOLVERLIST_H
#define	PASCSOLVERLIST_H

namespace pascinference {
	
	enum SolverType{ 
		SOLVER_AUTO, /* choose automatic solver */
		SOLVER_CG, /* Conjugate Gradient method (Unconstrained QP problem) */
		SOLVER_SPGQP /* Spectral Projected Gradient method of QP (QP problem with given projection) */
		
		};
	
}


#endif
