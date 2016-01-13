#include "gamma.h"

void Gamma::init(Data data, int K)
{
	int k;
	
	/* set input values */
	this->K = K;
	this->T = data.get_T();

	/* prepare array with gamma vectors */
	this->gamma_vecs = new GammaVector<Scalar>[this->K];

	/* alloc first vector */
	GammaVector<Scalar> D(this->T);
	/* set initial zero value to all vectors */
	D(all) = 0.0;
	for(k=0;k<this->K;k++){
		this->gamma_vecs[k] = D;
	}

}

void Gamma::finalize()
{
	delete []this->gamma_vecs;
}

void Gamma::prepare_random()
{
	int k,t;
	GammaVector<Scalar> gamma_sum(this->T);
		
	/* generate random data to gamma */
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: could be performed fully parallel
			this->gamma_vecs[k](t) = rand()/(double)(RAND_MAX);
		}
	}
	
	/* normalize gamma */
	/* at first sum the vectors */
	gamma_sum = this->gamma_vecs[0];
	for(k=1;k<this->K;k++){
		gamma_sum += this->gamma_vecs[k];
	}

	/* now divide the gamma by gamma_sum value */
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: could be performed fully parallel
			if(gamma_sum(t) == 0){
				/* maybe we generated only zeros */
				if(k == 0){
					this->gamma_vecs[k](t) = 1.0;
				} else {
					this->gamma_vecs[k](t) = 0.0;
				}	
			} else {
				this->gamma_vecs[k](t) = this->gamma_vecs[k](t)/gamma_sum(t);
			}
		}	
	}

}

void Gamma::prepare_uniform()
{
	int k;
	Scalar value;
	
	/* generate gamma = 1/K for all T */
	value = 1.0/(Scalar)this->K;
	for(k=0;k<this->K;k++){
		this->gamma_vecs[k](all) = value;
	}
}

//void Gamma::compute(QPSolver *qpsolver, Data data, Theta theta)
//{
	/* --- SOLVE OPTIMIZATION PROBLEM --- */
//	ierr = qpsolver->solve(); CHKERRQ(ierr);
//}
/*
void Gamma::compute_g(Vec g, Data *data, Theta *theta)
{
	void ierr;
	
	Vec g_part;
	PetscScalar *g_part_arr, *g_arr;
	int i,k;
		
	PetscFunctionBegin;

	ierr = VecDuplicate(this->gamma_vecs[0],&g_part); CHKERRQ(ierr);

	ierr = VecGetArray(g,&g_arr); CHKERRQ(ierr);
	for(k=0;k<this->dim;k++){
		ierr = this->compute_gk(g_part, data, theta, k); CHKERRQ(ierr);
		
		ierr = VecGetArray(g_part,&g_part_arr); CHKERRQ(ierr);
		for(i=0;i<this->local_size;i++){
			g_arr[k*this->local_size+i] = g_part_arr[i];
		}
		ierr = VecRestoreArray(g_part,&g_part_arr); CHKERRQ(ierr);
	}
	ierr = VecRestoreArray(g,&g_arr); CHKERRQ(ierr);

	ierr = VecAssemblyBegin(g); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(g); CHKERRQ(ierr);		

	ierr = VecDestroy(&g_part); CHKERRQ(ierr);
	
    PetscFunctionReturn(0); 
}
*/
//void Gamma::compute_gk(Vec g, Data *data, Theta *theta, int k)
//{
//	Vec *x_minus_Theta;
//	PetscScalar alpha,p;
	
//	PetscScalar *alphas;
	
//	int i;
//	int g_local_begin;
		
//	PetscFunctionBegin;

	/* get the ownership range of given g */
//	ierr = VecGetOwnershipRange(g, &g_local_begin, NULL); CHKERRQ(ierr);

	/* prepare array with coefficients */
//	ierr = PetscMalloc(data->get_dim()*sizeof(PetscScalar), &alphas); CHKERRQ(ierr);
//	for(i=0;i<data->get_dim();i++){
//		alphas[i] = 1.0;
//	}
	
	/* prepare array of vectors x_minus_Theta */
//	ierr = PetscMalloc(data->get_dim()*sizeof(Vec), &x_minus_Theta); CHKERRQ(ierr);
//	for(i=0;i<data->get_dim();i++){
//		ierr = VecDuplicate(data->data_vecs[i],&(x_minus_Theta[i])); CHKERRQ(ierr);
//	}

//	alpha = -1.0;
//	p = 2.0;
//	for(i=0;i<data->get_dim();i++){
		/* x_minus_Theta = Theta */
//		ierr = VecSet(x_minus_Theta[i],theta->theta_arr[k*data->get_dim()+i]); CHKERRQ(ierr);

		/* x_minus_Theta = x - Theta */
//		ierr = VecAYPX(x_minus_Theta[i], alpha, data->data_vecs[i]); CHKERRQ(ierr);

		/* x_minus_Theta = x_minus_Theta.^2 */
//		ierr = VecPow(x_minus_Theta[i], p); CHKERRQ(ierr);
//	}
	
	/* compute the sum of x_minus_Theta[:] */
//	ierr = VecSet(g,0.0); CHKERRQ(ierr);
//	ierr = VecMAXPY(g, data->get_dim(), alphas, x_minus_Theta); CHKERRQ(ierr);

	/* destroy temp vectors */
//	for(i=0;i<data->get_dim();i++){
//		ierr = VecDestroy(&(x_minus_Theta[i])); CHKERRQ(ierr);
//	}
//	ierr = PetscFree(alphas); CHKERRQ(ierr);
	
//    PetscFunctionReturn(0); 
//}

void Gamma::print()
{
	int k;
	std::ostringstream oss;
	std::ostringstream oss_values;
	
	Message_info("- gamma:");
	for(k=0;k<this->K;k++){
		oss << " - gamma[" << k << "] = ";
		oss_values << this->gamma_vecs[k];
		Message_info_values(oss.str(),oss_values.str());

		oss.str("");
		oss.clear();
		oss_values.str("");
		oss_values.clear();
	}
}

void Gamma::print(int nmb_of_spaces)
{
	int k,i; /* iterator */

	std::ostringstream oss_spaces;

	std::ostringstream oss;
	std::ostringstream oss_values;
	
	for(i=0;i<nmb_of_spaces;i++){
		oss_spaces << " ";
	}
	
	oss << oss_spaces.str() << "- gamma:";
	Message_info(oss.str());
	oss.str("");
	oss.clear();
	
	for(k=0;k<this->K;k++){
		oss << oss_spaces.str() << " - gamma[" << k << "] = ";
		oss_values << this->gamma_vecs[k];
		Message_info_values(oss.str(),oss_values.str());

		oss.str("");
		oss.clear();
		oss_values.str("");
		oss_values.clear();
	}
	
}

int Gamma::get_T()
{
	return this->T;
}

int Gamma::get_K()
{
	return this->K;
}
