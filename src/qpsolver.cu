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
	this->time_init = 0.0;
	this->time_stepsize = 0.0;
	this->time_fs = 0.0;
	this->time_total = 0.0;
	
}

/* prepare data which are constant */
void QPSolver::init(){
	/* the time for initialization is the part of total time, it is necessary to add it */
	timer.start(); 
	
	int T = this->get_T();
	int K = this->get_K();
	
	/* prepare RHS bs, gs, ds */

	/* alloc first vector */
	DataVector<Scalar> b(K*T);
	/* set initial zero value to all vectors */
	b(all) = 0.0;

	this->b = b;
	this->g = b;
	this->d = b;
	this->Ad = b;

	this->it_all = 0;
	this->hessmult_all = 0;

	this->time_total += timer.stop();
}

void QPSolver::finalize(){
	/* the time is the part of total time, it is necessary to add it */
	timer.start(); 

	/* clean the mess */

	this->time_total += timer.stop();
}

void QPSolver::compute_b(){

	this->gamma->compute_gk(&this->b, this->data, this->theta);
	this->b *= -1.0;

}

void QPSolver::solve(){
	timer.start(); /* add to time total in the end of solution */
	timer.start(); /* here starts the counter of time_init */

	/* algorithm parameters */
	int m = ALGORITHM_SPGQP_m;
	Scalar gamma = ALGORITHM_SPGQP_gamma;
	Scalar sigma2 = ALGORITHM_SPGQP_sigma2;
	Scalar eps = ALGORITHM_SPGQP_eps;
	int maxit = ALGORITHM_SPGQP_maxit;
	Scalar alphainit = 1.0/(this->eps_sqr*ALGORITHM_SPGQP_lambdaest);

	/* output performance */
	this->it = 0;
	this->hessmult = 0;

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
	get_projection(&(this->gamma->gamma_vec), this->get_K(), &this->time_projection);

	/* compute gradient, g = A*x-b */
	get_Ax_laplace(this->g,this->gamma->gamma_vec,&this->time_matmult); 
	this->hessmult += 1; /* there was muliplication by A */
	this->g -= this->b;
	
	/* compute function value */
	fx = this->get_function_value(this->gamma->gamma_vec, true);
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

	this->time_init = timer.stop(); /* here stop the initialization */
	
	/* main cycle */
	while(this->it < maxit){
		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		timer.start(); /* this is vector update */
		this->d = this->gamma->gamma_vec - alpha_bb*(this->g);
		this->time_update += timer.stop();

		/* d = P(d) */
		get_projection(&this->d, K, &this->time_projection);
		
		/* d = d - x */
		/* Ad = A*d */
		timer.start();
		this->d += -this->gamma->gamma_vec;
		this->time_update += timer.stop();

		get_Ax_laplace(this->Ad,this->d,&this->time_matmult);
		this->hessmult += 1; /* there was muliplication by A */

		/* dd = dot(d,d) */
		/* dAd = dot(Ad,d) */
		/* gd = dot(g,d) */
		dd = get_dot(this->d,this->d,&this->time_dot);
		dAd = get_dot(this->Ad,this->d,&this->time_dot);
		gd = get_dot(this->g,this->d,&this->time_dot);
		
		/* stopping criteria */
		if(dd < eps){
			break;
		}
		
		/* fx_max = max(fs) */
		timer.start(); /* manipulation with fs */
		fx_max = max(fs);	
		this->time_fs += timer.stop();
		
		/* compute step-size from A-condition */
		timer.start(); /* step-size timer */
		xi = (fx_max - fx)/dAd;
		beta_bar = -gd/dAd;
		beta_hat = gamma*beta_bar + sqrt(gamma*gamma*beta_bar*beta_bar + 2*xi);

		/* beta = min(sigma2,beta_hat) */
		if(beta_hat < sigma2){
			beta = beta_hat;
		} else {
			beta = sigma2;
		}
		this->time_stepsize += timer.stop();

		/* update approximation and gradient */
		/* x = x + beta*d */
		timer.start();/* this is vector update */
		this->gamma->gamma_vec += (this->d)*beta; 
		this->time_update += timer.stop();

		/* use recursive formula to compute gradient */
		/* g = g + beta*Ad */
		timer.start();
		this->g += (this->Ad)*beta;
		this->time_update += timer.stop();
		
		/* compute new function value using gradient */
		fx = this->get_function_value(this->gamma->gamma_vec, true);
		
		/* update fs */
		/* fs(1:end-1) = fs(2:end); */
		/* fs(end) = f;	*/
		timer.start(); /* manipulation with fs */
		if(m == 1){
			fs(0) = fx;
		} else {
			fs(0,m-2) = fs(1,m-1);
			fs(m-1) = fx;
		}
		this->time_fs += timer.stop();
		
		/* update BB step-size */
		timer.start(); /* step-size timer */
		alpha_bb = dd/dAd;
		this->time_stepsize += timer.stop();
		
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

	this->it_all += this->it;
	this->hessmult_all += this->hessmult;

	/* say goodbye */
	if(DEBUG_ALGORITHM_SAYHELLO){
		Message_info_main("\n- final info:");
		Message_info_time(" - time: \t\t",this->time_total);
		Message_info_value(" - it: \t\t\t",this->it);
		Message_info_value(" - hessmult: \t\t",this->hessmult);
		Message_info_value(" - final fx = \t\t",fx);
		Message_info("- SPGQP END ---------------------------------------------------------------");
	}

	/* very short info */
	if(DEBUG_ALGORITHM_BASIC){
		Message_info("  - SPGQP algorithm");
		Message_info_value("   - it    = ",this->it);
		Message_info_time("   - time  = ",this->time_total);

	}

	this->time_total += timer.stop();
}

Scalar QPSolver::get_function_value(){
	return this->get_function_value(this->gamma->gamma_vec);
}

Scalar QPSolver::get_function_value(GammaVector<Scalar> x){
	return this->get_function_value(x,false);
}

Scalar QPSolver::get_function_value(GammaVector<Scalar> x, bool use_gradient){
	Scalar fx = 0.0;

	if(use_gradient){
		/* use computed gradient in this->gs to compute function value */
		fx = 0.5*get_dot(this->g-this->b,x,&this->time_dot);
	} else {
		/* we have nothing - compute fx using full formula fx = 0.5*dot(A*x,x) - dot(b,x) */
		
		GammaVector<Scalar> Ax(this->get_T()*this->get_K());
		Scalar xAx, xb;

		get_Ax_laplace(Ax,x,&this->time_matmult);
		 
		xAx = get_dot(Ax,x,&this->time_dot);
		fx += 0.5*xAx;
		 
		xb = get_dot(x,this->b,&this->time_dot);
		fx -= xb;
	}	

	return fx;	
}


void QPSolver::print(){
	this->print(0);
}

void QPSolver::print(int nmb_of_spaces){
	int i,k;
	int K = this->get_K();
	int T = this->get_T();
	
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
		oss_values << this->b(k*T,(k+1)*T-1);
		Message_info_values(oss.str(),oss_values.str());	
		oss.str(""); oss.clear();
		oss_values.str(""); oss_values.clear();
	}

	oss << oss_spaces.str() << " - vector of unknowns x:";
	Message_info(oss.str());
	oss.str(""); oss.clear();
	for(k=0;k<K;k++){
		oss << oss_spaces.str() << "   x[" << k << "] = ";
		oss_values << this->gamma->gamma_vec(k*T,(k+1)*T-1);
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

int QPSolver::get_it_all(){
	return this->it_all;
}

int QPSolver::get_hessmult(){
	return this->hessmult;
}

int QPSolver::get_hessmult_all(){
	return this->hessmult_all;
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

double QPSolver::get_time_init(){
	return this->time_init;
}

double QPSolver::get_time_stepsize(){
	return this->time_stepsize;
}

double QPSolver::get_time_fs(){
	return this->time_fs;
}

double QPSolver::get_time_other(){
	return this->time_total - (this->time_projection + this->time_matmult + this->time_dot + this->time_update + this->time_init + this->time_stepsize + this->time_fs);
}

