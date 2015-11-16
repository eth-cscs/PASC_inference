#include "qpproblem.h"

PetscErrorCode QPproblem::init(PetscInt N, PetscInt N_local, PetscInt K){
	PetscErrorCode ierr; /* error handler */

	PetscFunctionBegin;

	/* set input values */
	this->N = N;
	this->N_local = N_local;
	this->K = K;

	/* initialize data for optimization problem */
	/* prepare hessian matrix */
	ierr = MatCreate(PETSC_COMM_WORLD,&this->A); CHKERRQ(ierr);
	ierr = MatSetSizes(this->A,this->K*this->N_local,this->K*this->N_local,this->K*this->N,this->K*this->N); CHKERRQ(ierr);
	ierr = MatSetFromOptions(this->A); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(this->A,3,NULL,3,NULL); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(this->A,3,NULL); CHKERRQ(ierr);
	ierr = MatSetOption(this->A,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);

	/* prepare RHS b */
	ierr = MatGetVecs(this->A,&this->b,NULL); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->b); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->b,"b"); CHKERRQ(ierr);		
	
	/* prepare matrix of equality constraints BE */
	ierr = MatCreate(PETSC_COMM_WORLD, &this->BE); CHKERRQ(ierr);
	ierr = MatSetSizes(this->BE,PETSC_DECIDE,this->K*this->N_local,this->N,this->K*this->N); CHKERRQ(ierr);
	ierr = MatSetFromOptions(this->BE); CHKERRQ(ierr);
	ierr = MatSetType(this->BE, MATAIJ); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(this->BE,1,NULL,this->K*this->N,NULL); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(this->BE,1,NULL); CHKERRQ(ierr);

	/* prepare equality matrix */
	ierr = MatCreate(PETSC_COMM_WORLD, &this->BE); CHKERRQ(ierr);
	ierr = MatSetSizes(this->BE,this->N_local, this->K*this->N_local, this->N, this->K*this->N); CHKERRQ(ierr);
	ierr = MatSetFromOptions(this->BE); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->BE,"BE"); CHKERRQ(ierr);
	ierr = MatSetType(this->BE, MATAIJ); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(this->BE,this->K,NULL,this->K,NULL); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(this->BE,this->K,NULL); CHKERRQ(ierr);	

	/* prepare equality vector */
	ierr = MatGetVecs(this->BE,NULL,&this->cE); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->cE,"cE"); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->cE); CHKERRQ(ierr);

	/* prepare vector of bound constraints lb */
	ierr = VecDuplicate(this->b,&this->lb); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->lb); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->lb,"lower bound"); CHKERRQ(ierr);

	/* prepare solution vector x */
	ierr = VecDuplicate(this->b,&this->x); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->x); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->x,"x"); CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::finalize(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	/* clean the mess */
	ierr = MatDestroy(&this->A); CHKERRQ(ierr);
	ierr = VecDestroy(&this->b); CHKERRQ(ierr);
	ierr = MatDestroy(&this->BE); CHKERRQ(ierr);
	ierr = VecDestroy(&this->cE); CHKERRQ(ierr);
	ierr = VecDestroy(&this->lb); CHKERRQ(ierr);
	ierr = VecDestroy(&this->x); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}

PetscErrorCode QPproblem::assemble_A(){
	PetscErrorCode ierr;
	PetscInt k,i;

	PetscFunctionBegin;

	/* fill hessian matrix */
	for(k=0;k<this->K;k++){
		for(i=0;i<this->N;i++){
			/* first row */
			if(i == 0){
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i + 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
			}
			/* common row */
			if(i > 0 && i < this->N-1){
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i + 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i, 2.0, INSERT_VALUES); CHKERRQ(ierr);
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i - 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
			}
			/* last row */
			if(i == this->N-1){
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i - 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
				ierr = MatSetValue(this->A, k*this->N + i, k*this->N + i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
			}
		}	
	}
	/* Hessian matrix is filled and prepared */
	ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	
	/* A = 0.5*A (to obtain 1/2*x^T*A*x - b^T*x) */
	ierr = MatScale(this->A, 10*0.5); CHKERRQ(ierr); // TODO: eps_sqr
	
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::set_b(Vec b){
	PetscFunctionBegin;

	this->b = b;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::get_b(Vec *b){
	PetscFunctionBegin;

	*b = this->b;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::set_x(Vec x){
	PetscFunctionBegin;

	this->x = x;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::get_x(Vec *x){
	PetscFunctionBegin;

	*x = this->x;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::assemble_BE(){
	PetscErrorCode ierr;
	PetscInt k,i;

	PetscFunctionBegin;

	/* fill BE matrix */
	for(k=0;k<this->K;k++){
		/* fill eye(n),eye(n),eye(n) */
		for(i=0;i<this->N;i++){
			ierr = MatSetValue(this->BE, i, k*this->N + i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(this->BE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->BE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	{ /* Vaclav's feature */
/*		Mat BEt;
		TRY( FllopMatTranspose(BE,MAT_TRANSPOSE_EXPLICIT,&BEt) );
		TRY( MatDestroy(&BE) );
		TRY( FllopMatTranspose(BEt,MAT_TRANSPOSE_IMPLICIT,&BE) );
*/	}	

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::assemble_cE(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = VecSet(this->cE,1.0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(this->cE); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->cE); CHKERRQ(ierr);	

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::assemble_lb(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = VecSet(this->lb,0.0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(this->lb); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->lb); CHKERRQ(ierr);	

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::print(PetscViewer v){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = PetscViewerASCIIPrintf(v,"- QP optimization problem:\n"); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- Hessian matrix A:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = MatView(this->A,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIPrintf(v,"- right hand-side vector b:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = VecView(this->b,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- matrix of equality constraints BE:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = MatView(this->BE,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIPrintf(v,"- vector of equality constraints cE:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = VecView(this->cE,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- vector of bound constraints lb:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = VecView(this->lb,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- vector of unknowns x:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = VecView(this->x,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);


	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::solve_permon(){
	PetscErrorCode ierr;

	QP qp; /* qp problem */
	QPS qps; /* qp solver */

	PetscScalar normb; /* norm of b */
	PetscScalar rtol, atol, dtol; /* algorithm settings */
	PetscInt maxit; /* max number of iterations */

	PetscFunctionBegin;

	/* start the fun with Permon */
    FllopInitialize(NULL,NULL,(char*)0);	

	/* create QP */
	ierr = QPCreate(PETSC_COMM_WORLD,&qp); CHKERRQ(ierr);
	/* add loaded data to QP */
	ierr = QPSetInitialVector(qp,this->x); CHKERRQ(ierr);
	ierr = QPSetRhs(qp,this->b); CHKERRQ(ierr);
	ierr = QPSetOperator(qp, this->A, QP_SYM_SYMMETRIC); CHKERRQ(ierr);
	ierr = QPSetEq(qp,this->BE,this->cE); CHKERRQ(ierr);
	ierr = QPSetBox(qp,this->lb,NULL); CHKERRQ(ierr);
	/* create the QP solver (QPS) */
	ierr = QPSCreate(PETSC_COMM_WORLD, &qps); CHKERRQ(ierr);
	/* insert the QP problem into the solver */
	ierr = QPSSetQP(qps, qp); CHKERRQ(ierr);
	/* set default QPS options */
	ierr = VecNorm(this->b,NORM_2,&normb); CHKERRQ(ierr);
	
	rtol  = 0.001/normb;
	atol  = PETSC_DEFAULT;
	dtol  = PETSC_DEFAULT;
	maxit = PETSC_DEFAULT;
	TRY( QPSSetTolerances(qps, rtol, atol, dtol, maxit) );
	/* set QPS options from the options database */
	ierr = QPSSetFromOptions(qps); CHKERRQ(ierr);

	/* --- SOLVE OPTIMIZATION PROBLEM --- */
	ierr = QPSSolve(qps); CHKERRQ(ierr);

	/* get the solution vector */
	ierr = QPGetSolutionVector(qp, &this->x); CHKERRQ(ierr);

	/* destroy QP, QPS and other Permon stuff */
	ierr = QPDestroy(&qp); CHKERRQ(ierr);
	ierr = QPSDestroy(&qps); CHKERRQ(ierr);

	/* it is not necessary to do more fun with Permon */
	ierr = FllopFinalize(); CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}
