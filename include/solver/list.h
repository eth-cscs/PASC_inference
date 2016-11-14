#ifndef PASCSOLVERLIST_H
#define	PASCSOLVERLIST_H

namespace pascinference {
namespace solver {
	
	enum SolverType { 
		SOLVER_AUTO, /* choose automatic solver */
		SOLVER_CG, /* Conjugate Gradient method (Unconstrained QP problem) */
		SOLVER_SPGQP, /* GPU-Spectral Projected Gradient method of QP (QP problem with given projection) */
		SOLVER_PERMON /* permon solver */
	};


}	
} /* end of namespace */


#endif
