#include "qpsolver.h"

/* SPGQP SETTINGS */
#define ALGORITHM_SPGQP_m 1
#define ALGORITHM_SPGQP_gamma 0.9
#define ALGORITHM_SPGQP_sigma2 1.0
#define ALGORITHM_SPGQP_eps 0.0001
#define ALGORITHM_SPGQP_maxit 10000
#define ALGORITHM_SPGQP_lambdaest 4.0
#define DEBUG_ALGORITHM_BASIC true /* basic information about used algorithm and parameters */
#define DEBUG_ALGORITHM_SAYHELLO false /* information about used algorithm and parameters */
#define DEBUG_ALGORITHM_PRINTF false /* print object function value in every iteration */
#define DEBUG_ALGORITHM_PRINTFS false /* print vector of object functions in every iteration */
#define DEBUG_ALGORITHM_PRINTCOEFF true /* print computed coefficients in every iteration */


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

	int T = this->get_T();
	int K = this->get_K();
	
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
	A_sub *= this->eps_sqr;
	this->A_sub = A_sub;

	/* prepare RHS bs, gs, ds */
	this->bs = new GammaVector<Scalar>[K];
	this->gs = new GammaVector<Scalar>[K];
	this->ds = new GammaVector<Scalar>[K];
	this->Ads = new GammaVector<Scalar>[K];
	/* alloc first vector */
	DataVector<Scalar> b(T);
	/* set initial zero value to all vectors */
	b(all) = 0.0;
	for(k=0;k<K;k++){
		this->bs[k] = b;
		this->gs[k] = b;
		this->ds[k] = b;
		this->Ads[k] = b;
	}

}

void QPSolver::finalize(){
	/* clean the mess */
	delete []this->bs;
	delete []this->gs;
	delete []this->ds;
	delete []this->Ads;
}

void QPSolver::compute_b(){
	int k;
	for(k=0;k<this->gamma->get_K();k++){ // TODO: parallel
		this->gamma->compute_gk(&(this->bs[k]), this->data, this->theta, k);
		this->bs[k] *= -1.0;
	}
}

void QPSolver::solve(){
	/* algorithm parameters */
	int m = ALGORITHM_SPGQP_m;
	Scalar gamma = ALGORITHM_SPGQP_gamma;
	Scalar sigma2 = ALGORITHM_SPGQP_sigma2;
	Scalar eps = ALGORITHM_SPGQP_eps;
	int maxit = ALGORITHM_SPGQP_maxit;
	Scalar alphainit = 1.0/(this->eps_sqr*ALGORITHM_SPGQP_lambdaest);

	/* output performance */
	int it = 0;
	int hess_mult = 0;
	Scalar comp_time;
	timer.start();

	int k; /* iterator through clusters */
	int K = this->get_K(); /* number of clusters */
	Scalar fx; /* function value */
	GammaVector<Scalar> fs(m); /* store function values for generalized A-condition */
	Scalar fx_max; /* max(fs) */
	Scalar xi, beta_bar, beta_hat,beta; /* for A-condition */
	Scalar dd; /* dot(d,d) */
	Scalar gd; /* dot(g,d) */
	Scalar dAd; /* dot(Ad,d) */
	Scalar alpha_bb; /* BB step-size */
	
	// TODO: temp
	Scalar fx_orig;
	Scalar fx_final;
	
	
	/* compute and set new RHS */
	/* b = -g(data,theta) */
	this->compute_b();

	fx_orig = this->get_function_value();
	
	/* project initial approximation to feasible set */
//	this->project(&(this->gamma->gamma_vecs));

	/* compute gradient, g = A*x-b */
	hess_mult += 1;
	for(k=0;k<K;k++){ // TODO: parallel
		this->gs[k] = this->A_sub*this->gamma->gamma_vecs[k]; 
		this->gs[k] -= this->bs[k];
	}
	
//	std::cout << this->A_sub << std::endl;
	
//	std::cout << "g:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->gs[k] << std::endl;
//	}

	/* compute function value */
	fx = this->get_function_value();
	fs(all) = fx;

//	std::cout << "fs:" << fs << std::endl;
	
	/* initial step-size */
	alpha_bb = alphainit;

//	std::cout << "alpha_bb:" << alphainit << std::endl;

	/* print basic informations about algorithm */
	if(DEBUG_ALGORITHM_SAYHELLO){
		Message_info("- SPGQP BEGIN -------------------------------------------------------------");
		Message_info_main("- parameters:");
		Message_info_value(" - m = \t\t\t",m);
		Message_info_value(" - gamma = \t\t",gamma);
		Message_info_value(" - sigma2 = \t\t",sigma2);
		Message_info_value(" - alpha_init = \t",alphainit);
		Message_info_value(" - init fx = \t\t",fx);
		Message_info_value(" - eps = \t\t",eps);
		Message_info_value(" - maxit = \t\t",maxit);
	}
	
	/* main cycle */
	while(it < maxit){
		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		for(k = 0; k < K;k++){ // TODO: parallel
			this->ds[k] = this->gamma->gamma_vecs[k] - alpha_bb*this->gs[k];
		}

//	std::cout << "d1:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->ds[k] << std::endl;
//	}
		
		/* d = P(d) */
		this->project(&(this->ds));

//	std::cout << "d2:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->ds[k] << std::endl;
//	}
		
		/* d = d - x */
		/* Ad = A*d */
		/* dd = dot(d,d) */
		/* dAd = dot(Ad,d) */
		/* gd = dot(g,d) */
		dd = 0.0;
		dAd = 0.0;
		gd = 0.0;
		hess_mult+=1;
		for(k = 0; k < K;k++){ // TODO: parallel
			this->ds[k] -= this->gamma->gamma_vecs[k];
			this->Ads[k] = this->A_sub*this->ds[k];

			dd += dot(this->ds[k],this->ds[k]);
			dAd += dot(this->Ads[k],this->ds[k]);
			gd += dot(this->gs[k],this->ds[k]);
		}

//	std::cout << "d:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->ds[k] << std::endl;
//	}
//	std::cout << "Ad:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->Ads[k] << std::endl;
//	}
//	std::cout << "dd:" << dd << std::endl;
//	std::cout << "dAd:" << dAd << std::endl;
//	std::cout << "gd:" << gd << std::endl;
		
		/* stopping criteria */
		if(dd < eps){
			break;
		}
		
		/* fx_max = max(fs) */
		fx_max = max(fs);	
		
		/* compute step-size from A-condition */
		xi = (fx_max - fx)/dAd;
		beta_bar = -gd/dAd;
		beta_hat = beta_bar; //gamma*beta_bar + sqrt(gamma*gamma*beta_bar*beta_bar + 2*xi);

//	std::cout << "fmax:" << fx_max << std::endl;
//	std::cout << "xi:" << xi << std::endl;
//	std::cout << "beta_bar:" << beta_bar << std::endl;
//	std::cout << "beta_hat:" << beta_hat << std::endl;

		/* beta = min(sigma2,beta_hat) */
		if(beta_hat < sigma2){
			beta = beta_hat;
		} else {
			beta = sigma2;
		}

//	std::cout << "beta:" << beta << std::endl;
		
		/* update approximation and gradient */
		/* x = x + beta*d */
		/* g = g + beta*Ad */

		std::cout << "fx_before:" << this->get_function_value() << std::endl;
		for(k = 0; k < K;k++){ // TODO: parallel
			this->gamma->gamma_vecs[k] += beta*this->ds[k]; 
//			this->gs[k] += beta*this->Ads[k];

			this->gs[k] = this->A_sub*this->gamma->gamma_vecs[k];
			this->gs[k] += -this->bs[k];

		}
		std::cout << "fx_after:" << this->get_function_value() << std::endl;

		

//	std::cout << "x:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->gamma->gamma_vecs[k] << std::endl;
//	}

//	std::cout << "g:" << std::endl;
//	for(k=0;k<K;k++){ // TODO: parallel
//		std::cout << this->gs[k] << std::endl;
//	}
		
		/* compute new function value */
		fx = this->get_function_value();

		
		/* update fs */
		/* fs(1:end-1) = fs(2:end); */
		/* fs(end) = f;	*/
		if(m == 1){
			fs(0) = fx;
		} else {
			fs(0,m-2) = fs(1,m-1);
			fs(m-1) = fx;
		}

//	std::cout << "fs:" << fs << std::endl;
		
		/* update BB step-size */
		alpha_bb = dd/dAd;

//	std::cout << "alpha_bb:" << alpha_bb << std::endl;
		
		/* print progress of algorithm */
		if(DEBUG_ALGORITHM_PRINTF || DEBUG_ALGORITHM_PRINTFS || DEBUG_ALGORITHM_PRINTCOEFF){
			std::cout << "\033[33mit = \033[0m" << it << std::endl;
		}

		if(DEBUG_ALGORITHM_PRINTF){
			std::cout << "\033[36m fx = \033[0m" << fx << std::endl;
		}

		if(DEBUG_ALGORITHM_PRINTFS){
			std::cout << "\033[36m fs = \033[0m" << fs << std::endl;
		}
		
		if(DEBUG_ALGORITHM_PRINTCOEFF){
			std::cout << "\033[36m dd = \033[0m" << dd << ",";
			std::cout << "\033[36m dAd = \033[0m" << dAd << ",";
			std::cout << "\033[36m gd = \033[0m" << gd << std::endl;
			
			std::cout << "\033[36m fx = \033[0m" << fx << ",";
			std::cout << "\033[36m fx_max = \033[0m" << fx_max << ",";
			std::cout << "\033[36m xi = \033[0m" << xi << std::endl;
			
			std::cout << "\033[36m beta_bar = \033[0m" << beta_bar << ",";
			std::cout << "\033[36m beta_hat = \033[0m" << beta_hat << ",";
			std::cout << "\033[36m beta = \033[0m" << beta << std::endl;
			
			std::cout << "\033[36m alpha_bb = \033[0m" << alpha_bb << std::endl;
			
		}
		
		/* increase iteration counter */
		it += 1;
	} /* main cycle end */

	comp_time = timer.stop();
	/* say goodbye */
	if(DEBUG_ALGORITHM_SAYHELLO){
		Message_info_main("\n- final info:");
		Message_info_time(" - time: \t\t",comp_time);
		Message_info_value(" - it: \t\t\t",it);
		Message_info_value(" - hess_mult: \t\t",hess_mult);
		Message_info_value(" - final fx = \t\t",fx);
		Message_info("- SPGQP END ---------------------------------------------------------------");
	}

	fx_final = this->get_function_value();

	/* very short info */
	if(DEBUG_ALGORITHM_BASIC){
		Message_info("  - SPGQP algorithm");
		Message_info_value("   - it    = ",it);
		Message_info_time("   - time  = ",comp_time);

		Message_info_value("   - f_orig   = ",fx_orig);
		Message_info_value("   - f_final  = ",fx_final);

	}

}

void QPSolver::project(GammaVector<Scalar> **x){
	int t,k;
	GammaVector<Scalar> x_sub(this->get_K());

	for(t=0;t<this->get_T();t++){	// TODO: this is the place, where the parallel impementation should make a point
		/* cut x_sub from x */
		for(k=0;k<this->get_K();k++){
			x_sub(k) = (*x)[k](t);
		}
		
		/* compute subprojection */
		this->project_sub(&x_sub);

		/* add x_sub back to x */
		for(k=0;k<this->get_K();k++){
			(*x)[k](t) = x_sub(k);
		}
	}

}

/* project x_sub to feasible set defined by equality and inequality constraints
 * sum(x_sub) = 1
 * x_sub >= 0
 */
void QPSolver::project_sub(GammaVector<Scalar> *x_sub){
	int n = this->get_K(); /* nmb of components of x_sub */
	int i;

	bool is_inside = true;
	
	/* control equality constraints */
	if(sum(*x_sub) != 1){ 
		is_inside = false;
	}
	
	/* control inequality constraints */
	for(i = 0; i < n; i++){ // TODO: could be performed parallely  
		if((*x_sub)(i) < 0.0){
			is_inside = false;
		}
	}

	/* if given point is not inside the feasible domain, then do projection */
	if(!is_inside){
		/* compute sorted x_sub */
		GammaVector<Scalar> y(n);
		for(i=0;i<n;i++){ // TODO: it is really necessary?
			y(i) = (*x_sub)(i); 
		}
		this->sort_bubble(&y);

		/* now perform analytical solution of projection problem */	
		Scalar t_hat = 0.0;
		i = n - 1;
		Scalar ti;

		while(i >= 1){
			ti = (sum(y(i,n-1)) - 1.0)/(Scalar)(n-i);
			if(ti >= y(i-1)){
				t_hat = ti;
				i = -1; /* break */
			} else {
				i = i - 1;
			}
		}

		if(i == 0){
			t_hat = (sum(y)-1.0)/(Scalar)n;
		}
    
		for(i = 0; i < n; i++){ // TODO: could be performed parallely  
			/* (*x_sub)(i) = max(*x_sub-t_hat,0); */
			ti = (*x_sub)(i) - t_hat;	
			if(ti > 0.0){
				(*x_sub)(i) = ti;
			} else {
				(*x_sub)(i) = 0.0;
			}
		}
	}
}

/* sort values of scalar vector */
void QPSolver::sort_bubble(GammaVector<Scalar> *x){
	int n = x->size();
	int i;
	int nnew;
	Scalar swap;

	while(n > 0){
		/* Iterate through x */
		nnew = 0;
		for(i=1;i<n;i++){
			/* Swap elements in wrong order */
			if ((*x)(i) < (*x)(i - 1)){
				swap = (*x)(i);
				(*x)(i) = (*x)(i-1);
				(*x)(i-1) = swap;
				nnew = i;
			}
        }
		n = nnew;
    }
}

Scalar QPSolver::get_function_value(){
	Scalar fx = 0.0;
	int k;
/*	for(k=0;k<this->get_K();k++){ // TODO: parallel
		 fx += 0.5*dot(this->gs[k]-this->bs[k],this->gamma->gamma_vecs[k]);
	}
*/ 

	GammaVector<Scalar> temp(this->get_T());

	for(k=0;k<this->get_K();k++){ // TODO: parallel
		 temp = this->A_sub*this->gamma->gamma_vecs[k];
		 fx += 0.5*dot(temp,this->gamma->gamma_vecs[k]);
		 fx -= dot(this->gamma->gamma_vecs[k],this->bs[k]);
	}	

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
	int K = this->get_K();
	
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
	Message_info_value(oss.str(),this->get_K());
	oss.str(""); oss.clear();
	
	oss << oss_spaces.str() << " - T = ";
	Message_info_value(oss.str(),this->get_T());
	oss.str(""); oss.clear();

	oss << oss_spaces.str() << " - dim = ";
	Message_info_value(oss.str(),this->get_dim());
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

int QPSolver::get_T(){
	return this->data->get_T();
}

int QPSolver::get_dim(){
	return this->data->get_dim();
}

int QPSolver::get_K(){
	return this->gamma->get_K();
}


