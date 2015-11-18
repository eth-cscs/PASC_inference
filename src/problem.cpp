#include "problem.h"

PetscRandom rnd_problem; /* random numbers generator */

PetscErrorCode get_problem_value1(PetscScalar *value1_out, PetscScalar *value2_out)
{
	PetscErrorCode ierr;
	PetscScalar covariance[4] = {0.001, 0.0, 0.0, 0.1};
	PetscScalar mu[2] = {0.25, 0};

	PetscFunctionBegin;

	ierr = my_mvnrnd_D2(mu, covariance, value1_out, value2_out); CHKERRQ(ierr);

//	*value1_out = 1;
//	*value2_out = 1;
	
    PetscFunctionReturn(0);  	
}

PetscErrorCode get_problem_value2(PetscScalar *value1_out, PetscScalar *value2_out)
{
	PetscErrorCode ierr;
	PetscScalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	PetscScalar mu[2] = {0.0, -0.5};

	PetscFunctionBegin;

	ierr = my_mvnrnd_D2(mu, covariance, value1_out, value2_out); CHKERRQ(ierr);
	
//	*value1_out = 2;
//	*value2_out = 2;
	
    PetscFunctionReturn(0);  	
}

PetscErrorCode get_problem_value3(PetscScalar *value1_out, PetscScalar *value2_out)
{
	PetscErrorCode ierr;
	PetscScalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	PetscScalar mu[2] = {0.0, 0.5};

	PetscFunctionBegin;

	ierr = my_mvnrnd_D2(mu, covariance, value1_out, value2_out); CHKERRQ(ierr);

//	*value1_out = 3;
//	*value2_out = 3;

    PetscFunctionReturn(0);  	
}

PetscErrorCode generate_problem(Data *data_out, PetscInt dataN)
{
	PetscErrorCode ierr;
	PetscInt j,j_global;
	Data data;
	PetscScalar *vec_array1,*vec_array2; /* for n=2 */
	PetscScalar random_value1, random_value2; 
	
	PetscFunctionBegin;

	/* prepare random generator */
	ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rnd_problem); CHKERRQ(ierr);
	ierr = PetscRandomSetType(rnd_problem,PETSCRAND); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rnd_problem); CHKERRQ(ierr);

	ierr = data.init(datan,dataN);

	/* generate random data */
	ierr = VecGetArray(data.data_vecs[0], &vec_array1); CHKERRQ(ierr); /* TODO: this is nasty, I know */
	ierr = VecGetArray(data.data_vecs[1], &vec_array2); CHKERRQ(ierr);

	for(j=0;j<data.get_local_size();j++){
		j_global = data.get_local_begin() + j;
		if(j_global >= 0 && j_global < data.get_global_size()/3){ 
			ierr = get_problem_value1(&random_value1, &random_value2);
		}
		if(j_global >= data.get_global_size()/3 && j_global < 2*data.get_global_size()/3){ 
			ierr = get_problem_value2(&random_value1, &random_value2);
		}
		if(j_global >= 2*data.get_global_size()/3 && j_global <= data.get_global_size()){ 
			ierr = get_problem_value3(&random_value1, &random_value2);
		}
		
		vec_array1[j] = random_value1;
		vec_array2[j] = random_value2;
	}

	ierr = VecRestoreArray(data.data_vecs[0], &vec_array1); CHKERRQ(ierr);
	ierr = VecRestoreArray(data.data_vecs[1], &vec_array2); CHKERRQ(ierr);

	ierr = VecAssemblyBegin(data.data_vecs[0]); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(data.data_vecs[0]); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(data.data_vecs[1]); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(data.data_vecs[1]); CHKERRQ(ierr);

	/* destroy the random generator */
	ierr = PetscRandomDestroy(&rnd_problem); CHKERRQ(ierr);

	*data_out = data;
    PetscFunctionReturn(0);  	
}

PetscErrorCode my_mvnrnd_D2(PetscScalar *mu, PetscScalar *covariance, PetscScalar *value1, PetscScalar *value2)
{
	PetscErrorCode ierr;

	PetscScalar L[4];
	PetscScalar r1, r2, r1n, r2n; 
	
	PetscReal R, c, s;

	PetscFunctionBegin;

    /* Compute normally distributed random numbers via Box-Muller */
    ierr = PetscRandomGetValueReal(rnd_problem, &r1);CHKERRQ(ierr);
    r1   = 1.-r1; /* to change from [0,1) to (0,1], which we need for the log */

    ierr = PetscRandomGetValueReal(rnd_problem, &r2);CHKERRQ(ierr);
    R = PetscSqrtReal(-2.*PetscLogReal(r1));
    c = PetscCosReal(2.*PETSC_PI*r2);
    s = PetscSinReal(2.*PETSC_PI*r2);

    /* compute normal distributed random values */
    r1n = R*c;
    r2n = R*s;
	
	/* choleski decomposition of SPD covariance matrix */
	L[0] = sqrt(covariance[0]);
	L[1] = 0;
	L[2] = covariance[1]/L[0];
	L[3] = sqrt(covariance[3] - L[2]*L[2]);

	/* compute output values */
	/* y = L*randn(2,1) + mean */
	*value1 = L[0]*r1n + L[1]*r2n + mu[0];
	*value2 = L[2]*r1n + L[3]*r2n + mu[1];

    PetscFunctionReturn(0);  	
	
}

