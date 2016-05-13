
namespace pascinference {

void Gamma::init(int dim, int T, int K)
{
	
	/* set input values */
	this->dim = dim;
	this->K = K;
	this->T = T;

	/* prepare gamma vector */
	GammaVector new_gamma_vec(this->T*this->K);
	new_gamma_vec(petscvector::all) = 0.0; // TODO: deal with all
	this->gamma_vec = new_gamma_vec;

	/* prepare QPSolver */
	QPSolver new_qpsolver;
	new_qpsolver.init(this->T, this->K, 10);
	if(DEBUG_MODE >= 10) new_qpsolver.print();

	this->qpsolver = new_qpsolver;

}

void Gamma::finalize()
{
	this->qpsolver.finalize();
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
	if(DEBUG_MODE >= 11) coutMaster << "  - gamma_sum_0 = " << gamma_sum << std::endl;

	for(k=1;k<this->K;k++){
		gamma_sum = gamma_sum + this->gamma_vec(k*T,(k+1)*T-1);

		if(DEBUG_MODE >= 11) coutMaster << "  - gamma_sum_" << k << " = " << gamma_sum << std::endl;
	}

	if(DEBUG_MODE >= 11) coutMaster << "  - gamma_sum = " << gamma_sum << std::endl;

	/* now divide the gamma by gamma_sum value */
	for(k=0;k<this->K;k++){
		for(t=0;t<this->T;t++){ // TODO: could be performed fully parallel
			if(gamma_sum(t) == 0.0){
				/* maybe we generated only zeros */
				if(k == 0){
					this->gamma_vec(k*T+t) = 1.0;
				} else {
					this->gamma_vec(k*T+t) = 0.0;
				}	
			} else {
				this->gamma_vec(k*T+t) /= gamma_sum(t);
			}
		}	
	}

	if(DEBUG_MODE >= 10) coutMaster << "  - generated gamma = " << this->gamma_vec << std::endl;

}

void Gamma::prepare_uniform()
{
	Scalar value;
	
	/* generate gamma = 1/K for all T */
	value = 1.0/(Scalar)this->K;
	this->gamma_vec(petscvector::all) = value; // TODO: deal with all

}

void Gamma::compute(DataVector data_vec, Theta theta)
{
	/* compute and set new RHS */
	Timer rhs_timer; // TODO: this timer should be elsewhere... btw. whole compute_gk function should be elsewhere
	rhs_timer.restart();
	rhs_timer.start();
	if(DEBUG_MODE >= 3) Message_info("   - preparing RHS");
		this->compute_gk(this->qpsolver.b, data_vec, theta);
	this->qpsolver.b *= -1.0;
	rhs_timer.stop();
	if(DEBUG_MODE >= 3) Message_info_time("   - RHS prepared in ", rhs_timer.get_value_sum());

	
	/* --- SOLVE OPTIMIZATION PROBLEM --- */
	this->qpsolver.solve(this->gamma_vec);
}

void Gamma::print()
{
	this->print(0);
}

void Gamma::print_timers()
{
	Message_info("  - gamma:");

	this->qpsolver.print_timers();

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

Scalar Gamma::get_function_value()
{
	return this->qpsolver.get_function_value(gamma_vec,false);
}

// TODO: temp
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

} /* end of namespace */
