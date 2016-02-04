#include "gamma.h"

void Gamma::init(int dim, int T, int K)
{
	
	/* set input values */
	this->dim = dim;
	this->K = K;
	this->T = T;

	/* prepare gamma vector */
	GammaVector new_gamma_vec(this->T*this->K);
	new_gamma_vec(all) = 0.0;
	this->gamma_vec = new_gamma_vec;

	/* prepare QPSolver */
	QPSolver new_qpsolver;
	new_qpsolver.init(this->gamma_vec, this->T, this->K, 100);
	if(DEBUG_MODE >= 3) new_qpsolver.print();

	this->qpsolver = new_qpsolver;

}

void Gamma::finalize()
{

}

void Gamma::prepare_random()
{
	int k,t;
	GammaVector gamma_sum(this->T);
		
	/* generate random data to gamma */
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: could be performed fully parallel
			this->gamma_vec(k*T+t) = rand()/(double)(RAND_MAX);
		}
	}
	
	/* normalize gamma */
	/* at first sum the vectors */
	gamma_sum = this->gamma_vec(0,T-1);
	for(k=1;k<this->K;k++){
		gamma_sum += this->gamma_vec(k*T,(k+1)*T-1);
	}

	/* now divide the gamma by gamma_sum value */
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: could be performed fully parallel
			if(gamma_sum(t) == 0){
				/* maybe we generated only zeros */
				if(k == 0){
					this->gamma_vec(k*T+t) = 1.0;
				} else {
					this->gamma_vec(k*T+t) = 0.0;
				}	
			} else {
				this->gamma_vec(k*T+t) = this->gamma_vec(k*T+t)/gamma_sum(t);
			}
		}	
	}

}

void Gamma::prepare_uniform()
{
	Scalar value;
	
	/* generate gamma = 1/K for all T */
	value = 1.0/(Scalar)this->K;
	this->gamma_vec(all) = value;

}

void Gamma::compute(DataVector data_vec, Theta theta)
{
	/* compute and set new RHS */
	this->compute_gk(this->qpsolver.b, data_vec, theta);
	this->qpsolver.b *= -1.0;
	
	/* --- SOLVE OPTIMIZATION PROBLEM --- */
	this->qpsolver.solve();
}

void Gamma::print()
{
	this->print(0);
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
	
	oss << oss_spaces.str() << "--- GAMMA ---";
	Message_info(oss.str());
	oss.str("");
	oss.clear();
	
	for(k=0;k<this->K;k++){
		oss << oss_spaces.str() << " - gamma[" << k << "] = ";
		oss_values << this->gamma_vec(k*T,(k+1)*T-1);
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

QPSolver Gamma::get_qpsolver()
{
	return this->qpsolver;
}

GammaVector Gamma::get_gamma_vec()
{
	return this->gamma_vec;
}


void Gamma::compute_gk(GammaVector& g, DataVector data_vec, Theta theta)
{
	int t,n,k;
	GammaVector temp(this->dim);
	
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: this could be performed parallely 
			for(n=0;n<this->dim;n++){
				temp(n) = data_vec(n*this->T+t) - theta.theta_vec(k*this->dim + n);
			} 
			g(k*this->T + t) = dot(temp,temp);
		}
	}
}

