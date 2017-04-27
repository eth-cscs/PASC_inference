/** @file taosolver.h
 *  @brief solve QP with TAO solvers
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_TAOSOLVER_H
#define	PASC_TAOSOLVER_H

#include <iostream>

#include "general/solver/qpsolver.h"
#include "general/data/qpdata.h"

#include "petsctao.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"

#define TAOSOLVER_DEFAULT_MAXIT 1000
#define TAOSOLVER_DEFAULT_EPS 1e-9
#define TAOSOLVER_USE_UPPERBOUND false
#define TAOSOLVER_DUMP false

namespace pascinference {
namespace solver {

/** \class TaoSolver
 *  \brief Interface with QP solvers implemented in TAO
 *
 *  For solving QP on closed convex set described by separable simplexes.
*/
template<class VectorBase>
class TaoSolver: public QPSolver<VectorBase> {
	private:
		/** \class TAOCtx
		 *  \brief The content of TAO solver for manipulation with object inside user-specific functions.
		 * 
		 */ 
		class TAOCtx {
			public:
				Mat A;	/**< Hessian matrix */
				Vec b;	/**< linear term */
				Vec x0;
				Mat B;	/**< matrix of equality constraints */
				Vec c;	/**< rhs of equality constraints */
				Vec cE;
				Vec temp_vec;	/**< auxiliary vector */
				Vec lb;
				Vec ub;
		};

		TAOCtx *userctx;
		static PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ctx);
		static PetscErrorCode FormHessian(Tao tao, Vec x, Mat H, Mat Hpre, void *ctx);
		static PetscErrorCode FormEqualityConstraints(Tao tao, Vec X, Vec CE, void* ctx);
		static PetscErrorCode FormEqualityJacobian(Tao tao, Vec X, Mat JE, Mat JEpre, void* ctx);

		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();

		Timer timer_solve; 			/**< total solution time of used algorithm */

		QPData<VectorBase> *qpdata; /**< data on which the solver operates */
		double gP; 					/**< norm of projected gradient */

		Tao taosolver;				/**< instance of TAO solver */

		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		/* temporary vectors used during the solution process */
		GeneralVector<VectorBase> *Ad; 		/**< A*d in cost function value computation */
//		Vec g_TAO; 		/**< gradient in internal TAO computation */
//		Vec eq_TAO; 	/**< value of equality constraints in internal TAO computation */

		/** @brief dump data of the solver
		 * 
		 * called before solve()
		 */
		void dump() const; 

		bool use_upperbound;		/**< use additional upper bound x<=1 */

	public:
		/** @brief general constructor
		* 
		*/
		TaoSolver();

		/** @brief constructor based on provided data of problem
		* 
		* @param new_qpdata data of quadratic program
		*/
		TaoSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		* 
		*/
		~TaoSolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printtimer(ConsoleOutput &output) const;
		void printshort(std::ostringstream &header, std::ostringstream &values) const;
		void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		std::string get_name() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

template<class VectorBase>
void TaoSolver<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("taosolver_maxit", &this->maxit, TAOSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("taosolver_eps", &this->eps, TAOSOLVER_DEFAULT_EPS);

	consoleArg.set_option_value("taosolver_use_upperbound", &this->use_upperbound, TAOSOLVER_USE_UPPERBOUND);	
	consoleArg.set_option_value("taosolver_dump", &this->dump_or_not, TAOSOLVER_DUMP);	
	
}


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
TaoSolver<VectorBase>::TaoSolver(){
	LOG_FUNC_BEGIN

	qpdata = NULL;
//	qp = NULL;
//	qps = NULL;
	
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;
	
	this->fx = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();
	
	userctx = new TAOCtx();
	
	LOG_FUNC_END
}

template<class VectorBase>
TaoSolver<VectorBase>::TaoSolver(QPData<VectorBase> &new_qpdata) {
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;	
	
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

	userctx = new TAOCtx();

	/* dissect QP objects from qpdata */
	// TODO: oh my god, this is not the best "general" way!
	BlockGraphSparseMatrix<PetscVector> *Abgs = dynamic_cast<BlockGraphSparseMatrix<PetscVector> *>(qpdata->get_A());
	userctx->A = Abgs->get_petscmatrix();
	double coeff = Abgs->get_coeff();
	TRYCXX( MatScale(userctx->A, coeff) );
	TRYCXX( MatAssemblyBegin(userctx->A,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(userctx->A,MAT_FINAL_ASSEMBLY) );	

	GeneralVector<PetscVector> *bg = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_b());
	userctx->b = bg->get_vector();

	GeneralVector<PetscVector> *x0g = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x0());
	userctx->x0 = x0g->get_vector();

	SimplexFeasibleSet_LinEqBound<PetscVector> *fs_lineqbound = dynamic_cast<SimplexFeasibleSet_LinEqBound<PetscVector> *>(qpdata->get_feasibleset());
	userctx->B = fs_lineqbound->get_B();
	userctx->c = fs_lineqbound->get_c();
	userctx->lb = fs_lineqbound->get_lb();

//	Mat BT;
//	TRYCXX( MatDuplicate(B, MAT_DO_NOT_COPY_VALUES, &BT) );
//	TRYCXX( MatTranspose(B, MAT_REUSE_MATRIX, &BT) );
//	TRYCXX( MatCreateTranspose(B, &BT) );

	/* prepare TAO solver */

	/* set type of solver:
		-#define TAOLMVM     "lmvm"		Limited Memory Variable Metric method for unconstrained minimization
		-#define TAONLS      "nls"		Newton's method with linesearch for unconstrained minimization
		-#define TAONTR      "ntr"		Newton's method with linesearch for unconstrained minimization
		-#define TAONTL      "ntl"	
		-#define TAOCG       "cg"		Newton's method with linesearch for unconstrained minimization
		*#define TAOTRON     "tron"		Newton Trust Region method for bound constrained minimization
		-#define TAOOWLQN    "owlqn"	Orthant-wise limited memory quasi-newton algorithm
		-#define TAOBMRM     "bmrm"		Bundle method for regularized risk minimization
		*#define TAOBLMVM    "blmvm"		Limited memory variable metric method for bound constrained minimization
		*#define TAOBQPIP    "bqpip"		Bounded quadratic interior point algorithm for quadratic optimization with box constraints
		*#define TAOGPCG     "gpcg"		Newton Trust Region method for quadratic bound constrained minimization
		*#define TAONM       "nm"		Gradient projected conjugate gradient algorithm is an active-set conjugate-gradient based method for bound-constrained minimization
		-#define TAOPOUNDERS "pounders"	Model-based algorithm pounder extended for nonlinear least squares
		-#define TAOLCL      "lcl"		Linearly constrained lagrangian method for pde-constrained optimization
		-#define TAOSSILS    "ssils"		Semi-smooth infeasible linesearch algorithm for solving complementarity constraints
		-#define TAOSSFLS    "ssfls"		Semi-smooth feasible linesearch algorithm for solving complementarity constraints
		-#define TAOASILS    "asils"		Active-set infeasible linesearch algorithm for solving complementarity constraints
		-#define TAOASFLS    "asfls"		Active-set feasible linesearch algorithm for solving complementarity constraints
		*#define TAOIPM      "ipm"		Interior point algorithm for generally constrained optimization. (Notes: This algorithm is more of a place-holder for future constrained optimization algorithms and should not yet be used for large problems or production code.)
	*/


	
//	TRYCXX( QPCreate(PETSC_COMM_WORLD, &qp) );
//	TRYCXX( QPSetOperator(qp, A) ); /* set stiffness matrix */
//	TRYCXX( QPSetRhs(qp, b) ); /* set righ hand-side vector */
//	TRYCXX( QPSetInitialVector(qp, x0) ); /* set initial approximation */

//	TRYCXX( QPAddEq(qp,B,c) ); /* add equality constraints Bx=c */
//	TRYCXX( QPSetBox(qp, lb, ub) ); /* add lowerbound */

	/* prepare permon QPS */
//	TRYCXX( QPSCreate(PETSC_COMM_WORLD, &qps) );
//	TRYCXX( QPSSetQP(qps, qp) ); /* Insert the QP problem into the solver. */
//	TRYCXX( QPSSetTolerances(qps, this->eps, this->eps, 1e12, this->maxit) ); /* Set QPS options from settings */
//	TRYCXX( QPSMonitorSet(qps,QPSMonitorDefault,NULL,0) ); /* Set the QPS monitor */
//	TRYCXX( QPTFromOptions(qp) ); /* Perform QP transforms */

//	TRYCXX( QPSSMALXESetRhoInitial(qps,1e5, QPS_ARG_DIRECT) );

//	TRYCXX( QPSSetFromOptions(qps) ); /* Set QPS options from the options database (overriding the defaults). */
//	TRYCXX( QPSSetUp(qps) ); /* Set up QP and QPS. */
	
	/* print some infos about QPS */
//	TRYCXX( QPSView(qps, PETSC_VIEWER_STDOUT_WORLD) );

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
TaoSolver<VectorBase>::~TaoSolver(){
	LOG_FUNC_BEGIN

	TRYCXX( TaoDestroy(&taosolver) );

	/* free temp vectors */
	free_temp_vectors();
	
	LOG_FUNC_END
}

template<class VectorBase>
PetscErrorCode TaoSolver<VectorBase>::FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ctx) {
	TAOCtx *userctx2 = (TAOCtx*)ctx;

	/* compute Ax - b */
	
	/* G = A*x */
	TRYCXX( MatMult(userctx2->A,X,G) );

	PetscScalar fx_quadratic, fx_linear;
	TRYCXX( VecDot(G,X, &fx_quadratic) );
	TRYCXX( VecDot(userctx2->b,X,&fx_linear) );
	*f = 0.5*fx_quadratic - fx_linear;

	/* G = -1*b + G */
	TRYCXX( VecAXPY(G,-1.0,userctx2->b) );

	/* fx = 0.5*x^T*A*x - b^T*x = 0.5*x^T*(g-b) */
//	*f = 10;
	
	//TODO: tests (temp)
	coutMaster << "TEST fx: " << *f << std::endl;
	PetscScalar normb;
	TRYCXX( VecNorm(userctx2->b,NORM_2, &normb) );
	coutMaster << "TEST normb: " << normb << std::endl;

	PetscScalar normc;
	TRYCXX( VecNorm(userctx2->c,NORM_2, &normc) );
	coutMaster << "TEST normc: " << normc << std::endl;

	PetscScalar normcE;
	TRYCXX( VecNorm(userctx2->cE,NORM_2, &normcE) );
//	TRYCXX( VecNorm(CE,NORM_2, &normcE) );
	coutMaster << "TEST normcE: " << normcE << std::endl;

	PetscScalar normlb;
	TRYCXX( VecNorm(userctx2->lb,NORM_2, &normlb) );
	coutMaster << "TEST normlb: " << normlb << std::endl;

//	PetscScalar normub;
//	TRYCXX( VecNorm(userctx2->ub,NORM_2, &normub) );
//	coutMaster << "TEST normub: " << normub << std::endl;


	Vec v;
	TRYCXX( VecDuplicate(userctx2->c,&v) );
	TRYCXX( MatMult(userctx2->B,X,v) );
	TRYCXX( VecAXPY(v,-1.0,userctx2->c) );
	double mynorm;
	TRYCXX( VecNorm(v, NORM_2, &mynorm) );
	coutMaster << "++++ norm(B*x - c) = " << mynorm << std::endl;
	TRYCXX( VecDestroy(&v));
	
//	TRYCXX( VecScale(G, -1.0) );
	
	//TRYCXX( VecView(G, PETSC_VIEWER_STDOUT_WORLD) );

	return(0);
}

template<class VectorBase>
PetscErrorCode TaoSolver<VectorBase>::FormHessian(Tao tao, Vec x, Mat H, Mat Hpre, void *ctx) {
//	TAOCtx *userctx2 = (TAOCtx*)ctx;

//	H = userctx2->A;
//	Hpre = userctx2->A;
//	TRYCXX( MatView(H, PETSC_VIEWER_STDOUT_WORLD) );

	return(0);
}

template<class VectorBase>
PetscErrorCode TaoSolver<VectorBase>::FormEqualityConstraints(Tao tao, Vec X, Vec CE, void* ctx) {
	TAOCtx *userctx2 = (TAOCtx*)ctx;

	/* CE = B*x - c */

	/* CE = B*x */
	TRYCXX( MatMult(userctx2->B,X,CE) );

	/* CE = -1*c + CE */
	TRYCXX( VecAXPY(CE,-1.0,userctx2->c) );

	return(0);
}

template<class VectorBase>
PetscErrorCode TaoSolver<VectorBase>::FormEqualityJacobian(Tao tao, Vec X, Mat JE, Mat JEpre, void* ctx) {
//	TAOCtx *userctx2 = (TAOCtx*)ctx;

//	JE = userctx2->B;
//	JEpre = userctx2->B;

	return(0);
}


/* prepare temp_vectors */
template<class VectorBase>
void TaoSolver<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<VectorBase> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */
	Ad = new GeneralVector<VectorBase>(*pattern);	
	
//	Vec b = qpdata->get_b()->get_vector();
//	TRYCXX( VecDuplicate( b, &g_TAO) );
	SimplexFeasibleSet_LinEqBound<PetscVector> *fs_lineqbound = dynamic_cast<SimplexFeasibleSet_LinEqBound<PetscVector> *>(qpdata->get_feasibleset());
	userctx->c = fs_lineqbound->get_c();
	TRYCXX( VecDuplicate( userctx->c, &(userctx->cE)) ); //TODO: destroy in free temp vectors
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void TaoSolver<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(Ad);
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void TaoSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
/*	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	output <<  " - m:          " << m << std::endl;
	output <<  " - gamma:      " << gamma << std::endl;
	output <<  " - sigma1:     " << sigma1 << std::endl;
	output <<  " - sigma2:     " << sigma2 << std::endl;
	output <<  " - alphainit:  " << alphainit << std::endl;
*/	
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}

	output.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
/*	output_local <<  " - maxit:      " << this->maxit << std::endl;
	output_local <<  " - eps:        " << this->eps << std::endl;
	output_local <<  " - debugmode: " << this->debugmode << std::endl;

	output_local <<  " - m:          " << m << std::endl;
	output_local <<  " - gamma:      " << gamma << std::endl;
	output_local <<  " - sigma1:     " << sigma1 << std::endl;
	output_local <<  " - sigma2:     " << sigma2 << std::endl;
	output_local <<  " - alphainit:  " << alphainit << std::endl;

	output_local.synchronize();
*/
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output_global, output_local);
		coutMaster.pop();
	}

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  " - it: " << std::setw(6) << this->it_last << ", ";
	output <<  "hess mult: " << std::setw(6) << this->hessmult_last << ", ";
	output <<  "fx: " << std::setw(10) << this->fx << ", ";	
//	output <<  "norm(gP): " << std::setw(10) << this->gP << ", ";
	output <<  "used memory: " << std::setw(6) << MemoryCheck::get_virtual() << "%" << std::endl;

	output << " - ";
//	output <<  "t_solve = " << std::setw(10) << this->timer_solve.get_value_last() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - fx:           " << std::setw(25) << this->fx << std::endl;
	output <<  "      - fx_control:   " << std::setw(25) << fx_quadratic+fx_linear << ", log: " << std::setw(25) << log(fx_quadratic+fx_linear)/log(10) << std::endl;
	output <<  "      - fx_linear:    " << std::setw(25) << fx_linear << ", log: " << std::setw(25) << log(fx_linear)/log(10) << std::endl;
	output <<  "      - fx_quadratic: " << std::setw(25) << fx_quadratic << ", log: " << std::setw(25) << log(fx_quadratic)/log(10) << std::endl;
//	output <<  "      - norm(gP):     " << std::setw(25) << this->gP << ", log: " << std::setw(25) << log(this->gP)/log(10) << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it all =        " << this->it_sum << std::endl;
	output <<  " - hessmult all =  " << this->hessmult_sum << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =      " << this->timer_solve.get_value_sum() << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);
	std::streamsize ss = std::cout.precision();

	values << std::setprecision(17);

	header << "TAO it, ";
	values << this->it_last << ", ";

	header << "TAO hessmult, ";
	values << this->hessmult_last << ", ";

	header << "TAO t all, ";
	values << this->timer_solve.get_value_last() << ", ";

	header << "TAO fx, ";
	values << this->fx << ", ";

	header << "TAO fx_linear, ";
	values << fx_linear << ", ";

	header << "TAO fx_quadratic, ";
	values << fx_quadratic << ", ";

	values << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void TaoSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	header << "TAO_sum it, ";
	values << this->it_sum << ", ";

	header << "TAO_sum hessmult, ";
	values << this->hessmult_sum << ", ";

	header << "TAO_sum t all, ";
	values << this->timer_solve.get_value_sum() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
std::string TaoSolver<VectorBase>::get_name() const {
	std::string return_value = "TaoSolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

/* solve the problem */
template<class VectorBase>
void TaoSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); /* stop this timer in the end of solution */

	TRYCXX( TaoCreate(PETSC_COMM_WORLD, &taosolver) );
	TRYCXX( TaoSetType(taosolver, TAOIPM) ); //TODO: is this the only one TAOsolver type?
	TRYCXX( TaoSetInitialVector(taosolver, userctx->x0) ); /* set initial approximation */
	if(this->use_upperbound){
		TRYCXX( VecDuplicate(userctx->lb,&(userctx->ub)) );
		TRYCXX( VecSet(userctx->ub,1.0) ); //TODO: destroy?
		TRYCXX( TaoSetVariableBounds(taosolver, userctx->lb, userctx->ub) );
	} else {
		userctx->ub = PETSC_NULL;
		TRYCXX( TaoSetVariableBounds(taosolver, userctx->lb, PETSC_NULL) );
	}
	TRYCXX( TaoSetObjectiveAndGradientRoutine(taosolver, FormFunctionGradient, (void*)userctx) );
	TRYCXX( TaoSetHessianRoutine(taosolver,userctx->A,userctx->A,FormHessian,(void*)userctx) );
	TRYCXX( TaoSetEqualityConstraintsRoutine(taosolver,userctx->cE,FormEqualityConstraints,(void*)userctx) );
	TRYCXX( TaoSetJacobianEqualityRoutine(taosolver, userctx->B, userctx->B, FormEqualityJacobian,(void*)userctx) );
//	TRYCXX( TaoSetTolerances(taosolver, this->eps, this->eps, this->eps) );
	TRYCXX( TaoSetTolerances(taosolver, 0, 0, 0) );
	TRYCXX( TaoSetFromOptions(taosolver) );


	int it = -1;
	int hessmult = -1;
	double fx = std::numeric_limits<double>::max();

	/* dump data */
	if(this->dump_or_not){
		this->dump();
	}

	/* final preparation of the solver */
	TRYCXX( TaoSetFromOptions(taosolver) );

	/* set KSP */
	KSP ksp;
	PC pc;
	TRYCXX( TaoGetKSP(taosolver,&ksp) );
	TRYCXX( KSPGetPC(ksp,&pc) );

//	TRYCXX( PCSetType(pc,PCJACOBI) );
	TRYCXX( PCSetType(pc,PCLU) );
	TRYCXX( PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU) );

	TRYCXX( KSPSetType(ksp,KSPPREONLY) );
	TRYCXX( KSPSetFromOptions(ksp) );

	/* call TAO solver */
	TRYCXX( TaoSolve(taosolver) );
	
	/* get the solution vector */
	GeneralVector<PetscVector> *xg = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x());
	Vec x = xg->get_vector();
	TRYCXX( TaoGetSolutionVector(taosolver, &x) );

	TRYCXX( TaoGetIterationNumber(taosolver, &it) );

	/* control the equality constraint */
	SimplexFeasibleSet_LinEqBound<PetscVector> *fs_lineqbound = dynamic_cast<SimplexFeasibleSet_LinEqBound<PetscVector> *>(qpdata->get_feasibleset());
	Mat B = fs_lineqbound->get_B();
	Vec c = fs_lineqbound->get_c();
	Vec v;
	TRYCXX( VecDuplicate(c,&v) );
	TRYCXX( MatMult(B,x,v) );
	TRYCXX( VecAXPY(v,-1.0,c) );

	double mynorm;
	TRYCXX( VecNorm(v, NORM_2, &mynorm) );
	coutMaster << "++++ norm(B*x - c) = " << mynorm << std::endl;


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
template<class VectorBase>
void TaoSolver<VectorBase>::dump() const {
	LOG_FUNC_BEGIN

	coutMaster << "... dump data ..." << std::endl;

	//TODO: write this function!

	LOG_FUNC_END
}


}
} /* end namespace */

#endif
