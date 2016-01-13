#include "qpsolver.h"

/* constructor */
QPSolver::QPSolver(Data* data, Gamma *gamma, Theta *theta, Scalar eps_sqr){
	this->data = data;
	this->gamma = gamma;
	this->theta = theta;
	this->eps_sqr = eps_sqr;	
}

/* prepare data which are constant */
void QPSolver::init(){
	int t,k;

	int T = this->data->get_T();
	int K = this->gamma->get_K();
	
	/* prepare block of hessian matrix */
	GammaMatrix<Scalar> A_sub(T,T); 
	for(t=0;t<T;t++){
		/* first row */
		if(t == 0){
			A_sub(t,t) = 1.0;
			A_sub(t,t+1) = -1.0;
		}
		/* common row */
		if(t > 0 && t < T-1){
			A_sub(t,t+1) = -1.0;
			A_sub(t,t) = 2.0;
			A_sub(t,t-1) = -1.0;
		}
		/* last row */
		if(t == T-1){
			A_sub(t,t-1) = -1.0;
			A_sub(t,t) = 1.0;
		}
	}	
	A_sub *= 0.5*this->eps_sqr;
	this->A_sub = A_sub;

	/* prepare RHS bs, gs */
	this->bs = new GammaVector<Scalar>[K];
	this->gs = new GammaVector<Scalar>[K];
	/* alloc first vector */
	DataVector<Scalar> b(T);
	/* set initial zero value to all vectors */
	b(all) = 0.0;
	for(k=0;k<K;k++){
		this->bs[k] = b;
		this->gs[k] = b;
	}

}

void QPSolver::finalize(){
	delete []this->bs;
	delete []this->gs;
}

void QPSolver::solve(){

}

Scalar QPSolver::get_function_value(){
	Scalar fx = 0.0;
	
	return fx;
}

void QPSolver::print(){
	this->print(0);
}

void QPSolver::print(int nmb_of_spaces){
	this->print(nmb_of_spaces, true);
}

void QPSolver::print(int nmb_of_spaces, bool print_A_sub){
	int i,k;
	int K = this->gamma->get_K();
	
	std::ostringstream oss_spaces;

	std::ostringstream oss;
	std::ostringstream oss_values;
	
	for(i=0;i<nmb_of_spaces;i++){
		oss_spaces << " ";
	}
	
	oss << oss_spaces.str() << "- QP optimization problem";
	Message_info(oss.str());
	oss.str(""); oss.clear();

	oss << oss_spaces.str() << " - K = ";
	Message_info_value(oss.str(),this->gamma->get_K());
	oss.str(""); oss.clear();
	
	oss << oss_spaces.str() << " - T = ";
	Message_info_value(oss.str(),this->data->get_T());
	oss.str(""); oss.clear();

	oss << oss_spaces.str() << " - dim = ";
	Message_info_value(oss.str(),this->data->get_dim());
	oss.str(""); oss.clear();

	if(print_A_sub){
		oss << oss_spaces.str() << " - block of Hessian matrix Asub:";
		Message_info(oss.str());
		oss.str(""); oss.clear();
		oss_values << this->A_sub;
		Message_info_values(oss.str(),oss_values.str());	
		oss_values.str(""); oss_values.clear();
	}

	oss << oss_spaces.str() << " - right hand-side vector b:";
	Message_info(oss.str());
	oss.str(""); oss.clear();
	for(k=0;k<K;k++){
		oss << oss_spaces.str() << "   b[" << k << "] = ";
		oss_values << this->bs[k];
		Message_info_values(oss.str(),oss_values.str());	
		oss.str(""); oss.clear();
		oss_values.str(""); oss_values.clear();
	}

	oss << oss_spaces.str() << " - vector of unknowns x:";
	Message_info(oss.str());
	oss.str(""); oss.clear();
	for(k=0;k<K;k++){
		oss << oss_spaces.str() << "   x[" << k << "] = ";
		oss_values << this->gamma->gamma_vecs[k];
		Message_info_values(oss.str(),oss_values.str());	
		oss.str(""); oss.clear();
		oss_values.str(""); oss_values.clear();
	}

}
