#include "qpproblemprojectionstep.h"

/* constructor */
QPproblemProjectionstep::QPproblemProjectionstep(Data *data, Gamma *gamma, Theta *theta, PetscScalar eps_sqr) : QPproblem(data, gamma,theta, eps_sqr) {
	this->N = this->gamma->get_global_size();
	this->N_local = this->gamma->get_local_size();
	this->K = this->gamma->get_dim();
	
	/* -0.99/lambda_max */
	this->stepsize = -0.99/(this->eps_sqr*4.0);

}

PetscErrorCode QPproblemProjectionstep::init(){
	PetscErrorCode ierr; /* error handler */
	PetscInt k;

	PetscFunctionBegin;

	/* initialize data for optimization problem */
	/* prepare block of hessian matrix */
	ierr = MatCreate(PETSC_COMM_WORLD,&this->Asub); CHKERRQ(ierr);
	ierr = MatSetSizes(this->Asub,this->N_local,this->N_local,this->N,this->N); CHKERRQ(ierr);
	ierr = MatSetFromOptions(this->Asub); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(this->Asub,3,NULL,3,NULL); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(this->Asub,3,NULL); CHKERRQ(ierr);
	ierr = MatSetOption(this->Asub,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);

	/* prepare RHS bs, gs */
	ierr = PetscMalloc(this->K*sizeof(Vec), &this->bs); CHKERRQ(ierr);
	ierr = PetscMalloc(this->K*sizeof(Vec), &this->gs); CHKERRQ(ierr);
	/* allocate the first vector, all other will be same */
	ierr = MatGetVecs(this->Asub,&this->bs[0],NULL); CHKERRQ(ierr);
	ierr = VecDuplicate(this->bs[0],&this->gs[0]); CHKERRQ(ierr);
	for(k=1;k<this->K;k++){
		ierr = VecDuplicate(this->bs[0],&this->bs[k]); CHKERRQ(ierr);
		ierr = VecDuplicate(this->gs[0],&this->gs[k]); CHKERRQ(ierr);
	}

	/* prepare temp vec */
	ierr = VecDuplicate(this->bs[0],&this->temp); CHKERRQ(ierr);
	ierr = VecDuplicate(this->bs[0],&this->temp2); CHKERRQ(ierr);

	/* we can immediately asseble some of QP objects which are independent of outer iterations */
	this->assemble_Asub();

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::finalize(){
	PetscErrorCode ierr;
	PetscInt i;

	PetscFunctionBegin;

	/* clean the mess */
	ierr = MatDestroy(&this->Asub); CHKERRQ(ierr);
	for(i=0;i<this->K;i++){
		ierr = VecDestroy(&this->bs[i]); CHKERRQ(ierr);
		ierr = VecDestroy(&this->gs[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(this->bs); CHKERRQ(ierr);
	ierr = PetscFree(this->gs); CHKERRQ(ierr);

	ierr = VecDestroy(&this->temp); CHKERRQ(ierr);
	ierr = VecDestroy(&this->temp2); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}

PetscErrorCode QPproblemProjectionstep::assemble_Asub(){
	PetscErrorCode ierr;
	PetscInt i;

	PetscFunctionBegin;

	/* fill the block of hessian matrix */
	for(i=0;i<this->N;i++){
		/* first row */
		if(i == 0){
			ierr = MatSetValue(this->Asub, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
			ierr = MatSetValue(this->Asub, i, i + 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
		}
		/* common row */
		if(i > 0 && i < this->N-1){
			ierr = MatSetValue(this->Asub, i, i + 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
			ierr = MatSetValue(this->Asub, i, i, 2.0, INSERT_VALUES); CHKERRQ(ierr);
			ierr = MatSetValue(this->Asub, i, i - 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
		}
		/* last row */
		if(i == this->N-1){
			ierr = MatSetValue(this->Asub, i, i - 1, -1.0, INSERT_VALUES); CHKERRQ(ierr);
			ierr = MatSetValue(this->Asub, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
		}
	}	
	/* Hessian matrix is filled and prepared */
	ierr = MatAssemblyBegin(this->Asub,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->Asub,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	
	/* A = 0.5*A (to obtain 1/2*x^T*A*x - b^T*x) */
	ierr = MatScale(this->Asub, this->eps_sqr*0.5); CHKERRQ(ierr);
	
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::set_bs(Vec *bs){
	PetscFunctionBegin;

	this->bs = bs;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::get_bs(Vec **bs){
	PetscFunctionBegin;

	*bs = this->bs;

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::print(PetscViewer v){
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = PetscViewerASCIIPrintf(v,"- QP optimization problem:\n"); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- block of Hessian matrix Asub:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	ierr = MatView(this->Asub,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIPrintf(v,"- right hand-side vector b:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
//	ierr = VecView(this->b,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

	ierr = PetscViewerASCIIPrintf(v,"- vector of unknowns x:\n"); CHKERRQ(ierr); 
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
//	ierr = VecView(this->x,v); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);


	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::solve(){
	PetscErrorCode ierr;
	PetscInt k;

	PetscFunctionBegin;

	/* compute gradient */
	ierr = this->compute_gradient(); CHKERRQ(ierr);

	/* --- PREPARE DATA FOR OPTIMIZATION PROBLEM --- */
	/* set new RHS, b = -g */
	for(k=0;k<this->K;k++){
		ierr = this->gamma->compute_gk(this->bs[k], this->data, this->theta, k); CHKERRQ(ierr); // TODO: move to model?
		ierr = VecScale(this->bs[k], -1.0); CHKERRQ(ierr);
		ierr = VecAssemblyBegin(this->bs[k]); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(this->bs[k]); CHKERRQ(ierr);	
	}

	// TODO: x = gamma

	/* --- SOLVE OPTIMIZATION PROBLEM --- */

	/* x = x - alpha*g */
	for(k=0;k<this->K;k++){
		ierr = VecAXPY(this->gamma->gamma_vecs[k], this->stepsize, this->gs[k]); CHKERRQ(ierr);
	}

	/* x = P(x) */
	ierr = this->project(); CHKERRQ(ierr);

	/* --- SET SOLUTION BACK TO GAMMA --- */
	// TODO: gamma = x

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::compute_gradient(){
	PetscErrorCode ierr;
	PetscInt k;

	PetscFunctionBegin;

	/* g = A*x - b */
	for(k=0;k<this->K;k++){
		ierr = MatMult(this->Asub, this->gamma->gamma_vecs[k], this->gs[k]); CHKERRQ(ierr);
		ierr = VecAXPY(this->gs[k], -1.0, this->bs[k]); CHKERRQ(ierr);
	}

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::get_function_value(PetscScalar *fx){
	PetscErrorCode ierr;
	PetscScalar value,my_sum;
	PetscInt k;
	
	PetscFunctionBegin;

	/* compute gradient */
	ierr = this->compute_gradient(); CHKERRQ(ierr);

	/* fx = 1/2*<g-b,x> */
	my_sum = 0.0;
	for(k=0;k<this->K;k++){
		ierr = VecCopy(this->gs[k],this->temp); CHKERRQ(ierr);
		ierr = VecAXPY(this->temp, -1.0, this->bs[k]); CHKERRQ(ierr);
		ierr = VecDot(this->gamma->gamma_vecs[k],this->temp,&value); CHKERRQ(ierr);
		my_sum += 0.5*value;
	}
	
	/* set return value */
	*fx = my_sum;

    PetscFunctionReturn(0);  
}


PetscErrorCode QPproblemProjectionstep::project(){
	PetscErrorCode ierr;
	PetscScalar norm_Bx;
	PetscScalar *alphas; /* = 1 */
	PetscScalar *betas; /* = -1/K */
	PetscInt k;
	PetscInt it; /* iteration counter */

	PetscFunctionBegin;

	/* prepare coeficients */
	ierr = PetscMalloc(this->K*sizeof(PetscScalar), &alphas); CHKERRQ(ierr);
	ierr = PetscMalloc(this->K*sizeof(PetscScalar), &betas); CHKERRQ(ierr);
	for(k=0;k<this->K;k++){
		alphas[k] = 1.0;
		betas[k] = -1.0/(PetscScalar)this->K;
	}	

	/* shift Bx=c to Bx=0 */
	/* gs[0] = -x_in */
	ierr = VecSet(this->gs[0],-1.0/(PetscScalar)this->K); CHKERRQ(ierr);
	for(k=0;k<this->K;k++){
		ierr = VecAXPY(this->gamma->gamma_vecs[k],1.0,this->gs[0]); CHKERRQ(ierr);
	}

	/* project to Bx=0 and x>=g */
	/* compute norm(Bx) as a stopping crit. */
	
	ierr = VecSet(this->temp2,0.0); CHKERRQ(ierr);
	ierr = VecMAXPY(this->temp2,this->K, alphas, this->gamma->gamma_vecs); CHKERRQ(ierr);
	ierr = VecNorm(this->temp2,NORM_2, &norm_Bx); CHKERRQ(ierr);
	
	it = 0;
	while(norm_Bx >= 0.0000001){ // TODO: this should be done in different way
		for(k=0;k<this->K;k++){
			/* project to Bx = 0 */
			ierr = VecCopy(this->gamma->gamma_vecs[k],this->temp); CHKERRQ(ierr);
			ierr = VecMAXPY(this->temp,this->K, betas, this->gamma->gamma_vecs); CHKERRQ(ierr);
			
			/* project to x>=g */
			ierr = VecPointwiseMax(this->gamma->gamma_vecs[k],this->temp,this->gs[0]); CHKERRQ(ierr);
		}
		
		/* compute norm(Bx) as a stopping crit. */
		ierr = VecSet(this->temp2,0.0); CHKERRQ(ierr);
		ierr = VecMAXPY(this->temp2,this->K, alphas, this->gamma->gamma_vecs); CHKERRQ(ierr);
		ierr = VecNorm(this->temp2,NORM_2, &norm_Bx); CHKERRQ(ierr);

		it += 1;
	}

	/* shift back Bx=0 to Bx=c */
	for(k=0;k<this->K;k++){
		ierr = VecAXPY(this->gamma->gamma_vecs[k],-1.0,this->gs[0]); CHKERRQ(ierr);
	}

	/* clean the mess */
	ierr = PetscFree(alphas); CHKERRQ(ierr);
	ierr = PetscFree(betas); CHKERRQ(ierr);

	/* print number of projection iterations; TODO: make it in different way */
	ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"- it_proj     = %d:\n",it); CHKERRQ(ierr);

    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblemProjectionstep::correct(PetscScalar increment){
	PetscFunctionBegin;
	
	this->stepsize = this->stepsize/2.0;
	
    PetscFunctionReturn(0);  
}
