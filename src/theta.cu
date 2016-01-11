#include "theta.h"

PetscErrorCode Theta::init(Data data, Gamma gamma){
	PetscErrorCode ierr;
	PetscInt i,k;
	
	PetscFunctionBegin;

	this->dim_data = data.get_dim();
	this->dim_gamma = gamma.get_dim();

	ierr = PetscMalloc(this->dim_data*this->dim_gamma*sizeof(PetscScalar),&this->theta_arr); CHKERRQ(ierr);
	/* set zeros to theta */
	for(k=0;k<this->dim_gamma;k++){
		for(i=0;i<this->dim_data;i++){
			this->theta_arr[k*this->dim_data + i] = 0.0;
		}
	}

    PetscFunctionReturn(0);  	
}

PetscErrorCode Theta::finalize()
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = PetscFree(this->theta_arr); CHKERRQ(ierr);

    PetscFunctionReturn(0);  	
}
	

PetscErrorCode Theta::compute(Data data, Gamma gamma){
	PetscErrorCode ierr;
	PetscScalar sum_gamma;
	PetscScalar gammaTx;
	
	PetscInt i,k;
		
	PetscFunctionBegin;

	for(k=0;k<gamma.get_dim();k++){

		/* compute sum of gamma[k] */
		ierr = VecSum(gamma.gamma_vecs[k],&sum_gamma); CHKERRQ(ierr);

		for(i=0;i<data.get_dim();i++){
			/* compute dot product */
			ierr = VecDot(gamma.gamma_vecs[k],data.data_vecs[i], &gammaTx); CHKERRQ(ierr);
			
			this->theta_arr[k*data.get_dim()+i] = gammaTx/sum_gamma;
		}
	}
	
    PetscFunctionReturn(0); 
}

PetscErrorCode Theta::print(PetscViewer v)
{
	PetscErrorCode ierr; /* error handler */
	PetscInt k,i; /* iterator */
	
	PetscFunctionBegin;
	
	ierr = PetscViewerASCIIPrintf(v,"- Theta:\n"); CHKERRQ(ierr);
	for(k=0;k<this->dim_gamma;k++){
		ierr = PetscViewerASCIIPrintf(v,"  - Theta_%d = ( ",k); CHKERRQ(ierr);
		for(i=0;i<this->dim_data;i++){
			ierr = PetscViewerASCIIPrintf(v,"%f ",this->theta_arr[k*this->dim_data+i]); CHKERRQ(ierr);
		}
		ierr = PetscViewerASCIIPrintf(v,")\n"); CHKERRQ(ierr);
	}

    PetscFunctionReturn(0);  		
}
