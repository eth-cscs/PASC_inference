#include "qpproblempermon.h"

/* constructor */
QPproblemPermon::QPproblemPermon(Data *data, Gamma *gamma, Theta *theta, PetscScalar eps_sqr) : QPproblem(data, gamma,theta, eps_sqr) {
	this->N = this->gamma->get_global_size();
	this->N_local = this->gamma->get_local_size();
	this->K = this->gamma->get_dim();
}

PetscErrorCode QPproblemPermon::init(){
	PetscErrorCode ierr; /* error handler */

	PetscFunctionBegin;

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

	/* prepare gradient vector g */
	ierr = VecDuplicate(this->b,&this->g); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->g); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->g,"g"); CHKERRQ(ierr);

	/* prepare temporary vector temp */
	ierr = VecDuplicate(this->b,&this->temp); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->temp); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->temp,"temp"); CHKERRQ(ierr);

	/* prepare temporary vector temp2 */
	ierr = VecDuplicate(this->cE,&this->temp2); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->temp2); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)this->temp2,"temp2"); CHKERRQ(ierr);

	/* prepare projector to ker BE */
	ierr = MatCreate(PETSC_COMM_WORLD,&this->PBE); CHKERRQ(ierr);
	ierr = MatSetSizes(this->PBE,this->K*this->N_local,this->K*this->N_local,this->K*this->N,this->K*this->N); CHKERRQ(ierr);
	ierr = MatSetFromOptions(this->PBE); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(this->PBE,this->N,NULL,this->K*this->N,NULL); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(this->PBE,this->K*this->N_local,NULL); CHKERRQ(ierr);
	ierr = MatSetOption(this->PBE,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);

	/* we can immediately asseble some of QP objects which are independent of outer iterations */
	this->assemble_A();
	this->assemble_BE();
	this->assemble_cE();
	this->assemble_lb();
	this->assemble_PBE();

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::finalize(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	/* clean the mess */
	ierr = MatDestroy(&this->A); CHKERRQ(ierr);
	ierr = VecDestroy(&this->b); CHKERRQ(ierr);
	ierr = MatDestroy(&this->BE); CHKERRQ(ierr);
	ierr = VecDestroy(&this->cE); CHKERRQ(ierr);
	ierr = VecDestroy(&this->lb); CHKERRQ(ierr);
	ierr = VecDestroy(&this->x); CHKERRQ(ierr);
	ierr = VecDestroy(&this->g); CHKERRQ(ierr);
	ierr = VecDestroy(&this->temp); CHKERRQ(ierr);
	ierr = MatDestroy(&this->PBE); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}

PetscErrorCode QPproblemPermon::assemble_A(){
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
	ierr = MatScale(this->A, this->eps_sqr*0.5); CHKERRQ(ierr);
	
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::assemble_PBE(){
	PetscErrorCode ierr;
	PetscInt k,i,j;

	PetscFunctionBegin;

	/* fill hessian matrix */
	for(k=0;k<this->K;k++){
		for(i=0;i<this->N;i++){
			ierr = MatSetValue(this->PBE, k*this->N + i, k*this->N + i, this->K-1.0, INSERT_VALUES); CHKERRQ(ierr);
			for(j=0;j<this->K;j++){
				if(j != k){
					ierr = MatSetValue(this->PBE, k*this->N + i, j*this->N + i, -1.0, INSERT_VALUES); CHKERRQ(ierr);
				}
			}
		}	
	}
	/* Hessian matrix is filled and prepared */
	ierr = MatAssemblyBegin(this->PBE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->PBE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	
	ierr = MatScale(this->PBE, 1.0/this->K); CHKERRQ(ierr);	
    
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::set_b(Vec b){
	PetscFunctionBegin;

	this->b = b;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::get_b(Vec *b){
	PetscFunctionBegin;

	*b = this->b;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::set_x(Vec x){
	PetscFunctionBegin;

	this->x = x;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::get_x(Vec *x){
	PetscFunctionBegin;

	*x = this->x;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::assemble_BE(){
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

PetscErrorCode QPproblemPermon::assemble_cE(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = VecSet(this->cE,1.0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(this->cE); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->cE); CHKERRQ(ierr);	

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::assemble_lb(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = VecSet(this->lb,0.0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(this->lb); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->lb); CHKERRQ(ierr);	

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::print(PetscViewer v){
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

	ierr = PetscViewerASCIIPrintf(v,"- projector to Ker BE:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = MatView(this->PBE,v); CHKERRQ(ierr);
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

PetscErrorCode QPproblemPermon::solve(){
	PetscErrorCode ierr;

	QP qp; /* qp problem */
	QPS qps; /* qp solver */

	PetscScalar normb; /* norm of b */
	PetscScalar rtol, atol, dtol; /* algorithm settings */
	PetscInt maxit; /* max number of iterations */

	PetscScalar *x_arr; /* array with local values of solution x */
	PetscScalar *gamma_arr;
	
	PetscInt k,i;
	PetscInt x_begin, x_end;

	PetscFunctionBegin;

	/* --- PREPARE DATA FOR OPTIMIZATION PROBLEM --- */

	/* set new RHS, b = -g */
	ierr = this->gamma->compute_g(this->b, this->data, this->theta); CHKERRQ(ierr);
	ierr = VecScale(this->b, -1.0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(this->b); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->b); CHKERRQ(ierr);	

	/* prepare initial vector from actual gamma, TODO: move this to Gamma */
	for(k=0;k<this->K;k++){
		ierr = VecGetArray(this->gamma->gamma_vecs[k],&gamma_arr); CHKERRQ(ierr);
		for(i=0;i<this->N_local;i++){
			ierr = VecSetValue(this->x, k*this->gamma->get_global_size()+this->gamma->get_local_begin()+i, gamma_arr[i], INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecRestoreArray(this->gamma->gamma_vecs[k],&gamma_arr); CHKERRQ(ierr);
	}
	ierr = VecAssemblyBegin(this->x); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->x); CHKERRQ(ierr);	



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


	/* --- SET SOLUTION BACK TO GAMMA --- */

	// TODO: move following to Gamma 

	/* get local array from solution x */
	ierr = VecGetOwnershipRange(this->x, &x_begin, &x_end); CHKERRQ(ierr);
	ierr = VecGetArray(this->x,&x_arr); CHKERRQ(ierr);
	for(i = x_begin;i < x_end; i++){
		k = (PetscInt)floor(i/(PetscScalar)this->N);
		ierr = VecSetValue(this->gamma->gamma_vecs[k], (i - k*this->N), x_arr[i-x_begin],INSERT_VALUES);
	}
	ierr = VecRestoreArray(this->x,&x_arr); CHKERRQ(ierr);

	/* the values of vectors are prepared for fun */
	for(k=0;k<this->K;k++){
		ierr = VecAssemblyBegin(this->gamma->gamma_vecs[k]); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(this->gamma->gamma_vecs[k]); CHKERRQ(ierr);	
	}


    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::compute_gradient(){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	/* g = A*x - b */
	ierr = MatMult(this->A, this->x, this->g); CHKERRQ(ierr);
	ierr = VecAXPY(this->g, -1.0, this->b); CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemPermon::get_function_value(PetscScalar *fx){
	PetscErrorCode ierr;
	PetscScalar value;
	
	PetscFunctionBegin;

	/* compute gradient */
	ierr = this->compute_gradient(); CHKERRQ(ierr);

	/* fx = 1/2*<g-b,x> */
	ierr = VecCopy(this->g,this->temp); CHKERRQ(ierr);
	ierr = VecAXPY(this->temp, -1.0, this->b); CHKERRQ(ierr);
	ierr = VecDot(this->x,this->temp,&value); CHKERRQ(ierr);
	value = 0.5*value;
	
	/* set return value */
	*fx = value;

    PetscFunctionReturn(0);  
}


PetscErrorCode QPproblemPermon::project(){
	PetscErrorCode ierr;
	PetscScalar norm_Bx;

	PetscFunctionBegin;

	/* shift Bx=c to Bx=0 */
	/* g = -x_in */
	ierr = VecSet(this->g,-1.0/(PetscScalar)this->K); CHKERRQ(ierr);
	ierr = VecAXPY(this->x,1.0,this->g); CHKERRQ(ierr);

	/* project to Bx=0 and x>=g */
	/* compute norm(Bx) as a stopping crit. */
	ierr = MatMult(this->BE,this->x,this->temp2); CHKERRQ(ierr);
	ierr = VecNorm(this->temp2,NORM_2, &norm_Bx); CHKERRQ(ierr);
	
	while(norm_Bx >= 0.00001){ // TODO: this should be done in different way
		/* project to Bx = 0 */
		ierr = MatMult(this->PBE,this->x,this->temp); CHKERRQ(ierr);
		
		/* project to x>=g */
		ierr = VecPointwiseMax(this->x,this->temp,this->g); CHKERRQ(ierr);
		
		/* compute norm(Bx) as a stopping crit. */
		ierr = MatMult(this->BE,this->x,this->temp2); CHKERRQ(ierr);
		ierr = VecNorm(this->temp2,NORM_2, &norm_Bx); CHKERRQ(ierr);
	}

	/* shift back Bx=0 to Bx=c */
	ierr = VecAXPY(this->x,-1.0,this->g); CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}
