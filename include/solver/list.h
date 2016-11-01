#ifndef PASCSOLVERLIST_H
#define	PASCSOLVERLIST_H

namespace pascinference {
namespace solver {
	
	enum SolverType{ 
		SOLVER_AUTO, /* choose automatic solver */
		SOLVER_QP, /* general QP solver (i.e. solve by SOLVER_AUTO in QPSolver) */
		SOLVER_CG, /* Conjugate Gradient method (Unconstrained QP problem) */
		SOLVER_SPGQP /* Spectral Projected Gradient method of QP (QP problem with given projection) */
		
	};


}	
} /* end of namespace */


#endif
