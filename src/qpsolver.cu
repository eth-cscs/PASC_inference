#include "qpsolver.h"

/* SPGQP SETTINGS */
#define ALGORITHM_SPGQP_m 30
#define ALGORITHM_SPGQP_gamma 0.9
#define ALGORITHM_SPGQP_sigma2 1.0
#define ALGORITHM_SPGQP_eps 0.0001
#define ALGORITHM_SPGQP_maxit 10000
#define ALGORITHM_SPGQP_lambdaest 4.0
#define DEBUG_ALGORITHM_BASIC false /* basic information about used algorithm and parameters */
#define DEBUG_ALGORITHM_SAYHELLO false /* information about used algorithm and parameters */
#define DEBUG_ALGORITHM_PRINTF false /* print object function value in every iteration */
#define DEBUG_ALGORITHM_PRINTFS false /* print vector of object functions in every iteration */
#define DEBUG_ALGORITHM_PRINTCOEFF false /* print computed coefficients in every iteration */


/* constructor */
QPSolver::QPSolver(Data* data, Gamma *gamma, Theta *theta, Scalar eps_sqr){
	this->data = data;
	this->gamma = gamma;
	this->theta = theta;
	this->eps_sqr = eps_sqr;
	
	this->time_projection = 0.0;
	this->time_matmult = 0.0;
	this->time_dot = 0.0;
	this->time_update = 0.0;
	this->time_total = 0.0;
	
}

/* prepare data which are constant */
void QPSolver::init(){
	/* the time for initialization is the part of total time, it is necessary to add it */
	timer.start(); 
	
	int k;

	int T = this->get_T();
	int K = this->get_K();
	
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

	this->time_total += timer.stop();
}

void QPSolver::finalize(){
	/* the time is the part of total time, it is necessary to add it */
	timer.start(); 

	/* clean the mess */
	delete []this->bs;
	delete []this->gs;
	delete []this->ds;
	delete []this->Ads;

	this->time_total += timer.stop();
}

void QPSolver::compute_b(){
	int k;
	for(k=0;k<this->gamma->get_K();k++){ // TODO: parallel
		this->gamma->compute_gk(&(this->bs[k]), this->data, this->theta, k);
		this->bs[k] *= -1.0;
	}
}

void QPSolver::solve(){
	timer.start(); /* add to time total in the end of solution */

	/* algorithm parameters */
	int m = ALGORITHM_SPGQP_m;
	Scalar gamma = ALGORITHM_SPGQP_gamma;
	Scalar sigma2 = ALGORITHM_SPGQP_sigma2;
	Scalar eps = ALGORITHM_SPGQP_eps;
	int maxit = ALGORITHM_SPGQP_maxit;
	Scalar alphainit = 1.0/(this->eps_sqr*ALGORITHM_SPGQP_lambdaest);

	/* output performance */
	this->it = 0;
	this->hess_mult = 0;

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
	
	/* compute and set new RHS */
	/* b = -g(data,theta) */
	this->compute_b();

	/* project initial approximation to feasible set */
	this->project(&(this->gamma->gamma_vecs));

	/* compute gradient, g = A*x-b */
	for(k=0;k<K;k++){ // TODO: parallel
		get_Ax_laplace(&(this->gs[k]),this->gamma->gamma_vecs[k],&this->time_matmult); 
		this->gs[k] -= this->bs[k];
	}
	
	/* compute function value */
	fx = this->get_function_value(this->gamma->gamma_vecs);
	fs(all) = fx;
	
	/* initial step-size */
	alpha_bb = alphainit;

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
	while(this->it < maxit){
		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		for(k = 0; k < K;k++){ // TODO: parallel
			timer.start();
			this->ds[k] = this->gamma->gamma_vecs[k] - alpha_bb*(this->gs[k]);
			this->time_update += timer.stop();
		}

		/* d = P(d) */
		this->project(&(this->ds));
		
		/* d = d - x */
		/* Ad = A*d */
		/* dd = dot(d,d) */
		/* dAd = dot(Ad,d) */
		/* gd = dot(g,d) */
		dd = 0.0;
		dAd = 0.0;
		gd = 0.0;
		for(k = 0; k < K;k++){ // TODO: parallel
			timer.start();
			this->ds[k] += -this->gamma->gamma_vecs[k];
			this->time_update += timer.stop();

			get_Ax_laplace(&(this->Ads[k]),this->ds[k],&this->time_matmult);

			dd += get_dot(this->ds[k],this->ds[k],&this->time_dot);
			dAd += get_dot(this->Ads[k],this->ds[k],&this->time_dot);
			gd += get_dot(this->gs[k],this->ds[k],&this->time_dot);
		}
		
		/* stopping criteria */
		if(dd < eps){
			break;
		}
		
		/* fx_max = max(fs) */
		fx_max = max(fs);	
		
		/* compute step-size from A-condition */
		xi = (fx_max - fx)/dAd;
		beta_bar = -gd/dAd;
		beta_hat = gamma*beta_bar + sqrt(gamma*gamma*beta_bar*beta_bar + 2*xi);

		/* beta = min(sigma2,beta_hat) */
		if(beta_hat < sigma2){
			beta = beta_hat;
		} else {
			beta = sigma2;
		}

		/* update approximation and gradient */
		/* x = x + beta*d */
		/* g = g + beta*Ad */
		for(k = 0; k < K;k++){ // TODO: parallel
			timer.start();
			this->gamma->gamma_vecs[k] += (this->ds[k])*beta; 
			this->time_update += timer.stop();

			/* use recursive formula to compute gradient */
			timer.start();
			this->gs[k] += (this->Ads[k])*beta;
			this->time_update += timer.stop();

//			this->gs[k] = this->A_sub*this->gamma->gamma_vecs[k];
//			this->gs[k] -= this->bs[k];

		}
		
		/* compute new function value using gradient */
		fx = this->get_function_value(this->gamma->gamma_vecs, true);
		
		/* update fs */
		/* fs(1:end-1) = fs(2:end); */
		/* fs(end) = f;	*/
		if(m == 1){
			fs(0) = fx;
		} else {
			fs(0,m-2) = fs(1,m-1);
			fs(m-1) = fx;
		}
		
		/* update BB step-size */
		alpha_bb = dd/dAd;
		
		/* print progress of algorithm */
		if(DEBUG_ALGORITHM_PRINTF || DEBUG_ALGORITHM_PRINTFS || DEBUG_ALGORITHM_PRINTCOEFF){
			std::cout << "\033[33mit = \033[0m" << this->it << std::endl;
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
		this->it += 1;
	} /* main cycle end */

	this->time_total += timer.stop();

	/* say goodbye */
	if(DEBUG_ALGORITHM_SAYHELLO){
		Message_info_main("\n- final info:");
		Message_info_time(" - time: \t\t",this->time_total);
		Message_info_value(" - it: \t\t\t",this->it);
		Message_info_value(" - hess_mult: \t\t",this->hess_mult);
		Message_info_value(" - final fx = \t\t",fx);
		Message_info("- SPGQP END ---------------------------------------------------------------");
	}

	/* very short info */
	if(DEBUG_ALGORITHM_BASIC){
		Message_info("  - SPGQP algorithm");
		Message_info_value("   - it    = ",this->it);
		Message_info_time("   - time  = ",this->time_total);

	}

}

void QPSolver::project(GammaVector<Scalar> **x){
	timer.start(); /* add to projection time */

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

	this->time_projection += timer.stop();
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

Scalar QPSolver::get_function_value(GammaVector<Scalar> *x){
	return this->get_function_value(x,false);
}

Scalar QPSolver::get_function_value(GammaVector<Scalar> *x, bool use_gradient){
	Scalar fx;
	int k;

	if(use_gradient){
		/* use computed gradient in this->gs to compute function value */
		for(k=0;k<this->get_K();k++){ // TODO: parallel
			fx += 0.5*get_dot(this->gs[k]-this->bs[k],this->gamma->gamma_vecs[k],&this->time_dot);
		}
	} else {
		/* we have nothing - compute fx using full formula fx = 0.5*dot(A*x,x) - dot(b,x) */
		
		GammaVector<Scalar> Ax(this->get_T());
		Scalar xAx, xb;

		fx = 0.0;
		for(k=0;k<this->get_K();k++){ // TODO: parallel
			get_Ax_laplace(&Ax,x[k],&this->time_matmult);
		 
			xAx = get_dot(Ax,x[k],&this->time_dot);
			fx += 0.5*xAx;
		 
			xb = get_dot(x[k],this->bs[k],&this->time_dot);
			fx -= xb;
		}

	}	

	return fx;	
}


void QPSolver::print(){
	this->print(0);
}

void QPSolver::print(int nmb_of_spaces){
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

int QPSolver::get_it(){
	return this->it;
}

int QPSolver::get_hessmult(){
	return this->hess_mult;
}

double QPSolver::get_time_projection(){
	return this->time_projection;
}

double QPSolver::get_time_matmult(){
	return this->time_matmult;
}

double QPSolver::get_time_dot(){
	return this->time_dot;
}

double QPSolver::get_time_update(){
	return this->time_update;
}

double QPSolver::get_time_total(){
	return this->time_total;
}

double QPSolver::get_time_other(){
	return this->time_total - (this->time_projection + this->time_matmult + this->time_dot + this->time_update);
}

