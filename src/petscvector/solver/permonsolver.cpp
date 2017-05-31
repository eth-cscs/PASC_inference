#include "external/petscvector/solver/permonsolver.h"

namespace pascinference {
namespace solver {

#ifdef USE_PERMON

template<>
PermonSolver<PetscVector>::PermonSolver(QPData<PetscVector> &new_qpdata){
	LOG_FUNC_BEGIN

	this->qpdata = &new_qpdata;
	
	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();
	externalcontent->qpdata = this->qpdata;
	externalcontent->qp = NULL;
	externalcontent->qps = NULL;
	
	/* allocate temp vectors */
	allocate_temp_vectors();
	
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;

	this->fx = std::numeric_limits<double>::max();
	this->gP = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	

	/* dissect QP objects from qpdata */
	// TODO: oh my god, this is not the best "general" way!
	BlockGraphSparseMatrix<PetscVector> *Abgs = dynamic_cast<BlockGraphSparseMatrix<PetscVector> *>(qpdata->get_A());

	Mat A = Abgs->get_externalcontent()->A_petsc;
	double coeff = Abgs->get_coeff();
	TRYCXX( MatScale(A, coeff) );
	TRYCXX( MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY) );	

	GeneralVector<PetscVector> *bg = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_b());
	Vec b = bg->get_vector();

	GeneralVector<PetscVector> *x0g = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x0());
	Vec x0 = x0g->get_vector();

	SimplexFeasibleSet_LinEqBound<PetscVector> *fs_lineqbound = dynamic_cast<SimplexFeasibleSet_LinEqBound<PetscVector> *>(qpdata->get_feasibleset());
	Mat B = fs_lineqbound->get_externalcontent()->B;
	Vec c = fs_lineqbound->get_externalcontent()->c;
	Vec lb = fs_lineqbound->get_externalcontent()->lb;

	/* prepare permon QP */
	TRYCXX( QPCreate(PETSC_COMM_WORLD, &(externalcontent->qp)) );
	TRYCXX( QPSetOperator(externalcontent->qp, A) ); /* set stiffness matrix */
	TRYCXX( QPSetInitialVector(externalcontent->qp, x0) ); /* set initial approximation */

	TRYCXX( QPAddEq(externalcontent->qp,B,c) ); /* add equality constraints Bx=c */
	if(this->use_upperbound){
		Vec ub;
		TRYCXX( VecDuplicate(lb,&ub) );
		TRYCXX( VecSet(ub,1.0) ); //TODO: destroy?
		TRYCXX( QPSetBox(externalcontent->qp, lb, ub) ); /* add box constraints */
	} else {
		TRYCXX( QPSetBox(externalcontent->qp, lb, PETSC_NULL) ); /* add lowerbound */
	}
	
	/* print some infos about QPS */
	TRYCXX( QPView(externalcontent->qp, PETSC_VIEWER_STDOUT_WORLD) );

	LOG_FUNC_END
}

/* solve the problem */
template<>
void PermonSolver<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	/* get QP and QPS */
	QP qp = externalcontent->qp;
	QPS qps = externalcontent->qps;

	/* prepare new transformations, vector b has been changed */
	GeneralVector<PetscVector> *bg = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_b());
	Vec b = bg->get_vector();

	BlockGraphSparseMatrix<PetscVector> *Abgs = dynamic_cast<BlockGraphSparseMatrix<PetscVector> *>(qpdata->get_A());
	double coeff = Abgs->get_coeff();

	/* prepare permon QPS */
	TRYCXX( QPSCreate(PETSC_COMM_WORLD, &qps) );
	TRYCXX( QPSSetQP(qps, qp) ); /* Insert the QP problem into the solver. */
	TRYCXX( QPSMonitorSet(qps,QPSMonitorDefault,NULL,0) ); /* Set the QPS monitor */
	TRYCXX( QPSetRhs(qp, b) ); /* set righ hand-side vector */
	TRYCXX( QPTFromOptions(qp) ); /* Perform QP transforms */
	TRYCXX( QPSSetFromOptions(qps) ); /* Set QPS options from the options database (overriding the defaults). */
	TRYCXX( QPSSetTolerances(qps, this->eps, this->eps, 1e12, this->maxit) ); /* Set QPS options from settings */

	/* provide max eigenvalue to PERMON */
	if(use_lambdamax){
		TRYCXX( QPSSMALXESetOperatorMaxEigenvalue(qps, 1.99*4.0*coeff) ); //TODO: only for Laplace!
	}

	TRYCXX( QPSSetUp(qps) ); /* Set up QP and QPS. */

	/* dump data */
	if(this->dump_or_not){
		externalcontent->dump();
	}
	
	this->timer_solve.start(); /* stop this timer in the end of solution */

	int it = -1;
	int hessmult = -1;
	double fx = std::numeric_limits<double>::max();

	/* call permon solver */
	TRYCXX( QPSSolve(qps) );
	
	/* get the solution vector */
	GeneralVector<PetscVector> *xg = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x());
	Vec x = xg->get_vector();
	TRYCXX( QPGetSolutionVector(qp, &x) );

	TRYCXX( QPSGetIterationNumber(qps, &it) );

	this->it_sum += it;
	this->hessmult_sum += hessmult;
	this->it_last = it;
	this->hessmult_last = hessmult;

	this->fx = fx;
	this->timer_solve.stop();

	/* write info to log file */
	LOG_IT(it)
	LOG_FX(fx)

	LOG_FUNC_END
}


/* dump data of the solver */
void PermonSolver<PetscVector>::ExternalContent::dump() const {
	LOG_FUNC_BEGIN

	coutMaster << "... dump data ..." << std::endl;
	
	/* prepare viewer to save to files */
	PetscViewer mviewer;

	/* A */
	Mat A;
	TRYCXX( QPGetOperator(qp, &A) );
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "A.bin" , FILE_MODE_WRITE, &mviewer) );
	TRYCXX( MatView(A, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	
	/* b */
	Vec b;
	TRYCXX( QPGetRhs(qp, &b) );
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "b.bin" , FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(b, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );

	/* B,c */
	Mat B;
	Vec c;
//	TRYCXX( QPGetEq(qp, &B, &c) );
	SimplexFeasibleSet_LinEqBound<PetscVector> *fs_lineqbound = dynamic_cast<SimplexFeasibleSet_LinEqBound<PetscVector> *>(qpdata->get_feasibleset());
	B = fs_lineqbound->get_externalcontent()->B;
	c = fs_lineqbound->get_externalcontent()->c;
	TRYCXX( VecView(c, PETSC_VIEWER_STDOUT_WORLD) );
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "Beq.bin" , FILE_MODE_WRITE, &mviewer) );
	TRYCXX( MatView(B, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "ceq.bin" , FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(c, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );

	/* x0 */
	Vec x0;
//	TRYCXX( QPGetInitialVector(qp, &x0) ); /* not working */
	GeneralVector<PetscVector> *x0g = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x0());
	x0 = x0g->get_vector();
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "x0.bin" , FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x0, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );


	LOG_FUNC_END
}

template<> 
PermonSolver<PetscVector>::ExternalContent * PermonSolver<PetscVector>::get_externalcontent() const {
	return this->externalcontent;	
}


#endif

}
}
