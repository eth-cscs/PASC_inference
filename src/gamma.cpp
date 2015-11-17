#include "gamma.h"

PetscErrorCode Gamma::init(Data data, PetscInt dim)
{
	PetscErrorCode ierr;
	PetscInt i;
	
	PetscFunctionBegin;

	this->dim = dim;
	this->global_size = data.get_global_size();
	this->local_size = data.get_local_size();

	/* get MPI variables */
    MPI_Comm_size(PETSC_COMM_WORLD,&this->proc_n);
    MPI_Comm_rank(PETSC_COMM_WORLD,&this->proc_id);

	ierr = PetscMalloc(this->dim*sizeof(Vec), &this->gamma_vecs); CHKERRQ(ierr);

	/* allocate the first vector, all other will be same */
	ierr = VecCreate(PETSC_COMM_WORLD,&this->gamma_vecs[0]); CHKERRQ(ierr);
	ierr = VecSetSizes(this->gamma_vecs[0],this->local_size,this->global_size); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->gamma_vecs[0]); CHKERRQ(ierr);
	for(i=1;i<this->dim;i++){
		ierr = VecDuplicate(this->gamma_vecs[0],&this->gamma_vecs[i]); CHKERRQ(ierr);
	}

	/* set ownership range */
	ierr = VecGetOwnershipRange(this->gamma_vecs[0], &this->local_begin, &this->local_end); CHKERRQ(ierr);

    PetscFunctionReturn(0);  	
}

PetscErrorCode Gamma::set_QPproblem(QPproblem*){
	PetscFunctionBegin;
	/* init QP problem */
	
    PetscFunctionReturn(0);  	
}


PetscErrorCode Gamma::finalize()
{
	PetscErrorCode ierr;
	PetscInt i;
	
	PetscFunctionBegin;

	/* destroy array with gamma vectors */
	for(i=1;i<this->dim;i++){
		ierr = VecDestroy(&this->gamma_vecs[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(this->gamma_vecs); CHKERRQ(ierr);

	/* finalize QP problem */
	this->qpproblem->finalize();

    PetscFunctionReturn(0);  	
}

PetscErrorCode Gamma::prepare_random()
{
	PetscErrorCode ierr;

	PetscInt k,i;
	PetscScalar *gamma_arr;
	PetscScalar *gamma_sum_arr;

	Vec gamma_sum;
	
	PetscFunctionBegin;

	/* prepare random generator */
	ierr = PetscRandomCreate(PETSC_COMM_WORLD,&this->rnd); CHKERRQ(ierr);
	ierr = PetscRandomSetType(this->rnd,PETSCRAND); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(this->rnd); CHKERRQ(ierr);

	/* generate random data to gamma */
	for(i=0;i<this->dim;i++){
		ierr = VecSetRandom(this->gamma_vecs[i], this->rnd); CHKERRQ(ierr);
	}
	
	/* normalize gamma */
	/* at first sum the vectors */
	ierr = VecDuplicate(this->gamma_vecs[0],&gamma_sum); CHKERRQ(ierr);
	ierr = VecCopy(this->gamma_vecs[0], gamma_sum); CHKERRQ(ierr);
	for(k=1;k<this->dim;k++){
		ierr = VecAXPY(gamma_sum,1.0,this->gamma_vecs[k]); CHKERRQ(ierr);
	}
	/* now divide the gamma by gamma_sum value */
	ierr = VecGetArray(gamma_sum, &gamma_sum_arr); CHKERRQ(ierr);
	for(k=0;k<this->dim;k++){
		ierr = VecGetArray(gamma_vecs[k], &gamma_arr); CHKERRQ(ierr);
		for(i=0;i<this->local_size;i++){
			if(gamma_sum_arr[i] == 0){
				/* maybe we generated only zeros */
				if(k == 0){
					gamma_arr[i] = 1.0;
				} else {
					gamma_arr[i] = 0.0;
				}	
			} else {
				gamma_arr[i] = gamma_arr[i]/gamma_sum_arr[i];
			}
		}	
		ierr = VecRestoreArray(this->gamma_vecs[k], &gamma_arr); CHKERRQ(ierr);
	}
	ierr = VecRestoreArray(gamma_sum, &gamma_sum_arr); CHKERRQ(ierr);

	/* destroy used sum vector */
	ierr = VecDestroy(&gamma_sum); CHKERRQ(ierr);
	
	/* destroy the random generator */
	ierr = PetscRandomDestroy(&this->rnd); CHKERRQ(ierr);

	/* assemble new values in vectors */
	for(k=0;k<this->dim;k++){
		ierr = VecAssemblyBegin(this->gamma_vecs[k]); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(this->gamma_vecs[k]); CHKERRQ(ierr);
	}	

    PetscFunctionReturn(0);  	
}

PetscErrorCode Gamma::prepare_uniform()
{
	PetscErrorCode ierr;

	PetscInt k;
	PetscScalar value;
	
	PetscFunctionBegin;

	/* generate gamma = 1/K for all T */
	value = 1.0/(PetscScalar)this->dim;
	for(k=0;k<this->dim;k++){
		ierr = VecSet(this->gamma_vecs[k], value); CHKERRQ(ierr);
	}

	for(k=0;k<this->dim;k++){
		ierr = VecAssemblyBegin(this->gamma_vecs[k]); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(this->gamma_vecs[k]); CHKERRQ(ierr);
	}	

    PetscFunctionReturn(0);  	
}

PetscErrorCode Gamma::prepare_fixed()
{
	PetscErrorCode ierr;

	PetscInt k;
	
	PetscFunctionBegin;

	// gamma0
	ierr = VecSet(this->gamma_vecs[0], 1.0); CHKERRQ(ierr);
	ierr = VecSetValue(this->gamma_vecs[0],0,0.0, INSERT_VALUES);
	ierr = VecSetValue(this->gamma_vecs[0],this->global_size-1,0.0, INSERT_VALUES);

	// gamma1
	ierr = VecSet(this->gamma_vecs[1], 0.0); CHKERRQ(ierr);
	ierr = VecSetValue(this->gamma_vecs[1],0,1.0, INSERT_VALUES);
	ierr = VecSetValue(this->gamma_vecs[1],this->global_size-1,0.0, INSERT_VALUES);

	// gamma2
	ierr = VecSet(this->gamma_vecs[2], 0.0); CHKERRQ(ierr);
	ierr = VecSetValue(this->gamma_vecs[2],0,0.0, INSERT_VALUES);
	ierr = VecSetValue(this->gamma_vecs[2],this->global_size-1,1.0, INSERT_VALUES);

	for(k=0;k<this->dim;k++){
		ierr = VecAssemblyBegin(this->gamma_vecs[k]); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(this->gamma_vecs[k]); CHKERRQ(ierr);
	}	

    PetscFunctionReturn(0);  	
}

PetscErrorCode Gamma::compute(Data data, Theta theta)
{
	PetscErrorCode ierr; /* error handler */
	
	PetscFunctionBegin;

	/* PRINT QP PROBLEM */
	if(PRINT_DATA){ /* print quadratic problem objects */
		ierr = this->qpproblem->print(PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}	

	/* --- SOLVE OPTIMIZATION PROBLEM --- */
	ierr = this->qpproblem->solve(); CHKERRQ(ierr);

    PetscFunctionReturn(0); 
}

PetscErrorCode Gamma::compute_g(Vec g, Data *data, Theta *theta)
{
	PetscErrorCode ierr;
	
	Vec *x_minus_Theta;
	Vec g_part;
	PetscScalar alpha,p;
	
	PetscScalar *alphas;
	PetscScalar *g_part_arr;
	
	PetscInt i,k;
	PetscInt g_local_begin;
		
	PetscFunctionBegin;

	/* get the ownership range of given g */
	ierr = VecGetOwnershipRange(g, &g_local_begin, NULL); CHKERRQ(ierr);

	/* prepare array with coefficients */
	ierr = PetscMalloc(data->get_dim()*sizeof(PetscScalar), &alphas); CHKERRQ(ierr);
	for(i=0;i<data->get_dim();i++){
		alphas[i] = 1.0;
	}
	
	/* prepare array of vectors x_minus_Theta */
	ierr = PetscMalloc(data->get_dim()*sizeof(Vec), &x_minus_Theta); CHKERRQ(ierr);
	for(i=0;i<data->get_dim();i++){
		ierr = VecDuplicate(data->data_vecs[i],&(x_minus_Theta[i])); CHKERRQ(ierr);
	}

	ierr = VecDuplicate(this->gamma_vecs[0],&g_part); CHKERRQ(ierr);

	alpha = -1.0;
	p = 2.0;
	for(k=0;k<this->get_dim();k++){
		for(i=0;i<data->get_dim();i++){
			/* x_minus_Theta = Theta */
			ierr = VecSet(x_minus_Theta[i],theta->theta_arr[k*data->get_dim()+i]); CHKERRQ(ierr);

			/* x_minus_Theta = x - Theta */
			ierr = VecAYPX(x_minus_Theta[i], alpha, data->data_vecs[i]); CHKERRQ(ierr);

			/* x_minus_Theta = x_minus_Theta.^2 */
			ierr = VecPow(x_minus_Theta[i], p); CHKERRQ(ierr);

		}
		/* compute the sum of x_minus_Theta[:] */
		ierr = VecSet(g_part,0.0); CHKERRQ(ierr);
		ierr = VecMAXPY(g_part, data->get_dim(), alphas, x_minus_Theta); CHKERRQ(ierr);

		ierr = VecAssemblyBegin(g_part); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(g_part); CHKERRQ(ierr);	

		/* view g_part */
		if(PETSC_FALSE){
			ierr = PetscPrintf(PETSC_COMM_WORLD,"- g_part_%d:\n",k); CHKERRQ(ierr);
			ierr = VecView(g_part,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
		}

		/* include g_part into vector g */
		ierr = VecGetArray(g_part, &g_part_arr); CHKERRQ(ierr);
		for(i=0;i<this->local_size;i++){
			ierr = VecSetValue(g, k*this->global_size+this->local_begin+i, g_part_arr[i], INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecRestoreArray(g_part, &g_part_arr); CHKERRQ(ierr);
	}

	ierr = VecAssemblyBegin(g); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(g); CHKERRQ(ierr);	

	/* destroy temp vectors */
	for(i=0;i<data->get_dim();i++){
		ierr = VecDestroy(&(x_minus_Theta[i])); CHKERRQ(ierr);
	}
	ierr = VecDestroy(&g_part); CHKERRQ(ierr);
	ierr = PetscFree(alphas); CHKERRQ(ierr);
	
    PetscFunctionReturn(0); 
}

PetscErrorCode Gamma::compute_gk(Vec g, Data *data, Theta *theta, PetscInt k)
{
	PetscErrorCode ierr;
	
	Vec *x_minus_Theta;
	PetscScalar alpha,p;
	
	PetscScalar *alphas;
	
	PetscInt i;
	PetscInt g_local_begin;
		
	PetscFunctionBegin;

	/* get the ownership range of given g */
	ierr = VecGetOwnershipRange(g, &g_local_begin, NULL); CHKERRQ(ierr);

	/* prepare array with coefficients */
	ierr = PetscMalloc(data->get_dim()*sizeof(PetscScalar), &alphas); CHKERRQ(ierr);
	for(i=0;i<data->get_dim();i++){
		alphas[i] = 1.0;
	}
	
	/* prepare array of vectors x_minus_Theta */
	ierr = PetscMalloc(data->get_dim()*sizeof(Vec), &x_minus_Theta); CHKERRQ(ierr);
	for(i=0;i<data->get_dim();i++){
		ierr = VecDuplicate(data->data_vecs[i],&(x_minus_Theta[i])); CHKERRQ(ierr);
	}

	alpha = -1.0;
	p = 2.0;
	for(i=0;i<data->get_dim();i++){
		/* x_minus_Theta = Theta */
		ierr = VecSet(x_minus_Theta[i],theta->theta_arr[k*data->get_dim()+i]); CHKERRQ(ierr);

		/* x_minus_Theta = x - Theta */
		ierr = VecAYPX(x_minus_Theta[i], alpha, data->data_vecs[i]); CHKERRQ(ierr);

		/* x_minus_Theta = x_minus_Theta.^2 */
		ierr = VecPow(x_minus_Theta[i], p); CHKERRQ(ierr);
	}
	
	/* compute the sum of x_minus_Theta[:] */
	ierr = VecSet(g,0.0); CHKERRQ(ierr);
	ierr = VecMAXPY(g, data->get_dim(), alphas, x_minus_Theta); CHKERRQ(ierr);

	/* destroy temp vectors */
	for(i=0;i<data->get_dim();i++){
		ierr = VecDestroy(&(x_minus_Theta[i])); CHKERRQ(ierr);
	}
	ierr = PetscFree(alphas); CHKERRQ(ierr);
	
    PetscFunctionReturn(0); 
}

PetscErrorCode Gamma::print(PetscViewer v)
{
	PetscErrorCode ierr; /* error handler */
	PetscInt i; /* iterator */
	
	PetscFunctionBegin;
	
	ierr = PetscViewerASCIIPrintf(v,"- gamma:\n"); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	for(i=0;i<this->dim;i++){
			ierr = PetscViewerASCIIPrintf(v,"- gamma_%d:\n",i); CHKERRQ(ierr);
			ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
			ierr = VecView(this->gamma_vecs[i],v); CHKERRQ(ierr);
			ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
	}
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}


PetscInt Gamma::get_local_size()
{
	return this->local_size;
}

PetscInt Gamma::get_global_size()
{
	return this->global_size;
}

PetscInt Gamma::get_local_begin()
{
	return this->local_begin;
}

PetscInt Gamma::get_local_end()
{
	return this->local_end;
}

PetscInt Gamma::get_dim()
{
	return this->dim;
}

PetscErrorCode Gamma::get_objectfunc_value(PetscScalar *value)
{
	PetscErrorCode ierr;	

	PetscFunctionBegin;

	/* compute function value from QP problem */
	ierr = this->qpproblem->get_function_value(value); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}





