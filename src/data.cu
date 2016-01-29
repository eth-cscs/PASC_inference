#include "data.h"

void Data::init(int dim, int T)
{

	/* set input values */
	this->dim = dim;
	this->T = T;
	
	/* prepare array with data vectors */
	DataVector<Scalar> new_data_vector(this->dim*this->T);
	new_data_vector(all) = 0.0;
	this->data_vec = new_data_vector;

}

void Data::finalize()
{

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
	
	Message_info("-- DATA ---");
	for(i=0;i<this->dim;i++){
		oss << "- data[" << i << "]: ";
		oss_values << this->data_vec(i*this->T,(i+1)*this->T-1);
		Message_info_values(oss.str(),oss_values.str());

		oss.str("");
		oss.clear();
		oss_values.str("");
		oss_values.clear();
	}
}

/**
 * generate random values in data vector 
 * 
**/
void Data::generate(){
	// TODO: here implement implicit random data generator, for example put random numbers everywhere
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

/* ------------ KMEANS ------------- */
void Data_kmeans::generate(){
	int t;
	int T = this->T;
	Scalar random_value1, random_value2; 
	
	/* generate random data */
	// TODO: this operation could be performed fully independently, no necessary to use for
	for(t=0;t<T;t++){ /* go through time serie */
		if(t >= 0 && t < T/3.0){ 
			this->get_problem_value1(&random_value1, &random_value2);
		}
		if(t >= T/3.0 && t < 2.0*T/3.0){ 
			this->get_problem_value2(&random_value1, &random_value2);
		}
		if(t >= 2.0*T/3.0 && t <= T){ 
			this->get_problem_value3(&random_value1, &random_value2);
		}
		
		this->data_vec(t) = random_value1;
		this->data_vec(T+t) = random_value2;
	}

}

void Data_kmeans::get_problem_value1(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.001, 0.0, 0.0, 0.1};
	Scalar mu[2] = {0.25, 0};

	this->my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void Data_kmeans::get_problem_value2(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	Scalar mu[2] = {0.0, -0.5};

	this->my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void Data_kmeans::get_problem_value3(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	Scalar mu[2] = {0.0, 0.5};

	this->my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void Data_kmeans::my_mvnrnd_D2(Scalar *mu, Scalar *covariance, Scalar *value1, Scalar *value2)
{
	Scalar L[4];
	Scalar r1, r2, r1n, r2n; 
	
	Scalar R, c, s;

    /* Compute normally distributed random numbers via Box-Muller */
	r1 = rand()/(double)(RAND_MAX);
    r1   = 1.-r1; /* to change from [0,1) to (0,1], which we need for the log */

    r2 = rand()/(double)(RAND_MAX);
    R = sqrt(-2.*log(r1));
    c = cos(2.*M_PI*r2);
    s = sin(2.*M_PI*r2);

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
}
