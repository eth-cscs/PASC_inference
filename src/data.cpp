#include "data.h"

PetscErrorCode Data::init(PetscInt dim, PetscInt global_size)
{
	PetscErrorCode ierr; /* error handler */
	PetscInt i; /* iterator */

	PetscFunctionBegin;

	/* get MPI variables */
    MPI_Comm_size(PETSC_COMM_WORLD,&this->proc_n);
    MPI_Comm_rank(PETSC_COMM_WORLD,&this->proc_id);

	/* set input values */
	this->dim = dim;
	this->global_size = global_size;
	
	/* prepare array with data vectors */
	ierr = PetscMalloc(dim*sizeof(Vec), &this->data_vecs); CHKERRQ(ierr);

	/* alloc first vector */
	ierr = VecCreate(PETSC_COMM_WORLD,&this->data_vecs[0]); CHKERRQ(ierr);
	ierr = VecSetSizes(this->data_vecs[0],PETSC_DECIDE,this->global_size); CHKERRQ(ierr);
	ierr = VecSetFromOptions(this->data_vecs[0]); CHKERRQ(ierr);
	/* all others vectors in data_vecs will be the same as data_vecs[0] */
	for(i=1;i<this->dim;i++){
		ierr = VecDuplicate(this->data_vecs[0],&this->data_vecs[i]); CHKERRQ(ierr);
	}

	/* set ownership range */
	ierr = VecGetLocalSize(this->data_vecs[0],&this->local_size); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(this->data_vecs[0], &this->local_begin, &this->local_end); CHKERRQ(ierr);

    PetscFunctionReturn(0);  	
}

PetscErrorCode Data::finalize()
{
	PetscErrorCode ierr;
	PetscInt i;
	
	PetscFunctionBegin;

	for(i=1;i<this->dim;i++){
		ierr = VecDestroy(&this->data_vecs[i]); CHKERRQ(ierr);
	}

	ierr = PetscFree(this->data_vecs); CHKERRQ(ierr);

    PetscFunctionReturn(0);  	
}


PetscInt Data::get_local_size()
{
	return this->local_size;
}

PetscInt Data::get_global_size()
{
	return this->global_size;
}

PetscInt Data::get_local_begin()
{
	return this->local_begin;
}

PetscInt Data::get_local_end()
{
	return this->local_end;
}

PetscInt Data::get_dim()
{
	return this->dim;
}

PetscErrorCode Data::print(PetscViewer v)
{
	PetscErrorCode ierr; /* error handler */
	PetscInt i; /* iterator */
	
	PetscFunctionBegin;
	
	ierr = PetscViewerASCIIPrintf(v,"- generated data:\n"); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
	for(i=0;i<this->dim;i++){
		ierr = PetscViewerASCIIPrintf(v,"- data[%d]:\n",i); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPushTab(v); CHKERRQ(ierr);
		ierr = VecView(this->data_vecs[i],v); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);
	}
	ierr = PetscViewerASCIIPopTab(v); CHKERRQ(ierr);

    PetscFunctionReturn(0);  		
}

PetscErrorCode Data::get_covtrace(PetscScalar *covtrace)
{
	PetscErrorCode ierr; /* error handler */
	PetscInt j; /* iterator */
	PetscScalar xTx, xTx_all;
	
	PetscFunctionBegin;

	xTx_all = 0;
	/* assemble new values in vetors */
	for(j=0;j<this->dim;j++){
		ierr = VecDot(this->data_vecs[j],this->data_vecs[j],&xTx); CHKERRQ(ierr);
		xTx_all += xTx;
	}	

	*covtrace = xTx_all;
    PetscFunctionReturn(0);  		
}

