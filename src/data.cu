#include "data.h"

void Data::init(int dim, int T)
{
	int i; /* iterator */

	/* set input values */
	this->dim = dim;
	this->T = T;
	
	/* prepare array with data vectors */
	this->data_vecs = new DataVector<Scalar>[this->dim];

	/* alloc first vector */
	DataVector<Scalar> D(this->T);
	/* set initial zero value to all vectors */
	D(all) = 0.0;
	for(i=0;i<this->dim;i++){
		this->data_vecs[i] = D;
	}

}

void Data::finalize()
{
	delete []this->data_vecs;
}


int Data::get_dim()
{
	return this->dim;
}

int Data::get_T()
{
	return this->T;
}

void Data::print()
{
	int i; /* iterator */
	std::ostringstream oss;
	std::ostringstream oss_values;
	
	Message_info("- generated data:");
	for(i=0;i<this->dim;i++){
		oss << " - data[" << i << "]: ";
		oss_values << this->data_vecs[i];
		Message_info_values(oss.str(),oss_values.str());

		oss.str("");
		oss.clear();
		oss_values.str("");
		oss_values.clear();
	}
}

//void Data::get_covtrace(PetscScalar *covtrace)
//{
//	PetscInt j;
//	PetscScalar xTx, xTx_all;
	
//	xTx_all = 0;
//	/* assemble new values in vetors */
//	for(j=0;j<this->dim;j++){
//		ierr = VecDot(this->data_vecs[j],this->data_vecs[j],&xTx); CHKERRQ(ierr);
//		xTx_all += xTx;
//	}	

//	*covtrace = xTx_all;
//}

