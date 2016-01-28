#include "theta.h"

void Theta::init(Data data, Gamma gamma){
	this->dim = data.get_dim();
	this->K = gamma.get_K();

	ThetaVector<Scalar> D(this->dim*this->K);
	this->theta_vec = D;
	this->theta_vec(all) = 0.0;

}

void Theta::finalize()
{

}
	

void Theta::compute(Data data, Gamma gamma){
	Scalar sum_gamma;
	Scalar gammaTx;
	int dim = data.get_dim();
	int T = gamma.get_T();
	int K = gamma.get_K();
	
	int i,k;
		
	for(k=0;k<K;k++){

		/* compute sum of gamma[k] */
		
		sum_gamma = sum(gamma.gamma_vec(k*T,(k+1)*T-1));

		for(i=0;i<dim;i++){
			/* compute dot product */
			gammaTx = dot(gamma.gamma_vec(k*T,(k+1)*T-1),data.data_vecs[i]);
			
			this->theta_vec(k*dim+i) = gammaTx/sum_gamma;
		}
	}
}

void Theta::print()
{
	this->print(0);
}

void Theta::print(int nmb_of_spaces)
{
	int k,i; /* iterator */
	int start_idx, end_idx;

	std::ostringstream oss_spaces;

	std::ostringstream oss;
	std::ostringstream oss_values;
	
	for(i=0;i<nmb_of_spaces;i++){
		oss_spaces << " ";
	}
	
	oss << oss_spaces.str() << "- Theta:";
	Message_info(oss.str());
	oss.str("");
	oss.clear();
	
	for(k=0;k<this->K;k++){
		start_idx = k*this->dim;
		end_idx = (k+1)*this->dim-1;
		
		oss << oss_spaces.str() << " - Theta_" << k << " = ";
		oss_values << this->theta_vec(range(start_idx,end_idx));
		Message_info_values(oss.str(),oss_values.str());

		oss.str("");
		oss.clear();
		oss_values.str("");
		oss_values.clear();
	}
	
}