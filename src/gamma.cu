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

void Gamma::compute(QPSolver *qpsolver, Data data, Theta theta)
{
	/* --- SOLVE OPTIMIZATION PROBLEM --- */
	qpsolver->solve();
}

void Gamma::compute_gk(GammaVector<Scalar> *g, Data *data, Theta *theta, int k)
{
	int t,n;
	GammaVector<Scalar> temp(data->get_dim());
	
	for(t=0;t<this->T;t++){ // TODO: this could be performed parallely 
		for(n=0;n<data->get_dim();n++){
			temp(n) = data->data_vecs[n](t) - theta->theta_vec(k*data->get_dim() + n);
		} 
		(*g)(t) = dot(temp,temp);
	}
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
