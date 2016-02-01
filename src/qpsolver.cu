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

/* prepare data which are constant */
void QPSolver::init(GammaVector<Scalar> x, int T, int K, Scalar eps_sqr){
	this->x = x;
	this->T = T;
	this->K = K;
	this->eps_sqr = eps_sqr;
	
	this->timer_projection.restart();
	this->timer_matmult.restart();
	this->timer_dot.restart();
	this->timer_update.restart();
	this->timer_stepsize.restart();
	this->timer_fs.restart();
	this->timer_total.restart();
	

	/* the time for initialization is the part of total time, it is necessary to add it */
	this->timer_total.start(); 
	
	/* prepare RHS bs, gs, ds */

	/* alloc first vector */
	DataVector<Scalar> b(this->K*this->T);
	/* set initial zero value to all vectors */
	b(all) = 0.0;

	this->b = b;
	this->g = b;
	this->d = b;
	this->Ad = b;

	this->it_all = 0;
	this->hessmult_all = 0;

	this->timer_total.stop();
}

void QPSolver::finalize(){
	/* the time is the part of total time, it is necessary to add it */
	this->timer_total.start(); 

	/* clean the mess */

	this->timer_total.stop();
}

void QPSolver::solve(){
	this->timer_total.start(); /* stop this timer in the end of solution */


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
		Message_info_value(" - eps = \t\t",eps);
		Message_info_value(" - maxit = \t\t",maxit);
	}

	/* project initial approximation to feasible set */
	this->timer_projection.start();
	 get_projection(this->x, this->get_K());
	this->timer_projection.stop();

	/* compute gradient, g = A*x-b */
	this->timer_matmult.start();
	 get_Ax_laplace(this->g,this->x,this->eps_sqr); 
 	 this->hessmult += 1; /* there was muliplication by A */
	this->timer_matmult.stop();
	this->g -= this->b;
	
	/* compute function value - it has its own timer */
	this->timer_fs.start();
 	 fx = this->get_function_value(this->x, true);
	 fs(all) = fx;
	this->timer_fs.stop();

	
	/* main cycle */
	while(this->it < maxit){
		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		this->timer_update.start(); /* this is vector update */
		 this->d = this->x - alpha_bb*(this->g);
		this->timer_update.stop();

		/* d = P(d) */
		this->timer_projection.start();
		 get_projection(this->d, K);
		this->timer_projection.stop();
		
		/* d = d - x */
		/* Ad = A*d */
		this->timer_update.start();
		 this->d += -this->x;
		this->timer_update.stop();

		this->timer_matmult.start();
		 get_Ax_laplace(this->Ad,this->d,this->eps_sqr);
		 this->hessmult += 1; /* there was muliplication by A */
		this->timer_matmult.stop();

		/* dd = dot(d,d) */
		/* dAd = dot(Ad,d) */
		/* gd = dot(g,d) */
		this->timer_dot.start();
		 dd = get_dot(this->d,this->d);
		 dAd = get_dot(this->Ad,this->d);
		 gd = get_dot(this->g,this->d);
		this->timer_dot.stop();
		
		/* stopping criteria */
		if(dd < eps){
			break;
		}
		
		/* fx_max = max(fs) */
		this->timer_fs.start(); /* manipulation with fs */
		 fx_max = max(fs);	
		this->timer_fs.stop();
		
		/* compute step-size from A-condition */
		this->timer_stepsize.start(); /* step-size timer */
		 xi = (fx_max - fx)/dAd;
		 beta_bar = -gd/dAd;
		 beta_hat = gamma*beta_bar + sqrt(gamma*gamma*beta_bar*beta_bar + 2*xi);

		 /* beta = min(sigma2,beta_hat) */
		 if(beta_hat < sigma2){
			beta = beta_hat;
		 } else {
			beta = sigma2;
		 }
		this->timer_stepsize.stop();

		/* update approximation and gradient */
		this->timer_update.start();/* this is vector update */
		 this->x += (this->d)*beta; /* x = x + beta*d */
		 this->g += (this->Ad)*beta; /* g = g + beta*Ad */
		this->timer_update.stop();
		
		/* compute new function value using gradient */
		this->timer_fs.start();
		 fx = this->get_function_value(this->x, true);
		
		 /* update fs */
		 /* fs(1:end-1) = fs(2:end); */
		 /* fs(end) = f;	*/
		 if(m == 1){
			fs(0) = fx;
		 } else {
			fs(0,m-2) = fs(1,m-1);
			fs(m-1) = fx;
		 }
		this->timer_fs.stop();
		
		/* update BB step-size */
		this->timer_stepsize.start(); /* step-size timer */
		 alpha_bb = dd/dAd;
		this->timer_stepsize.stop();
		
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

	this->timer_total.stop();

	/* say goodbye */
	if(DEBUG_ALGORITHM_SAYHELLO){
		Message_info_main("\n- final info:");
		Message_info_time(" - time: \t\t",this->timer_total.get_value_last());
		Message_info_value(" - it: \t\t\t",this->it);
		Message_info_value(" - hessmult: \t\t",this->hessmult);
		Message_info_value(" - final fx = \t\t",fx);
		Message_info("- SPGQP END ---------------------------------------------------------------");
	}

	/* very short info */
	if(DEBUG_ALGORITHM_BASIC){
		Message_info("  - SPGQP algorithm");
		Message_info_value("   - it    = ",this->it);
		Message_info_time("   - time  = ",this->timer_total.get_value_last());

	}

}

Scalar QPSolver::get_function_value(){
	return this->get_function_value(this->x);
}

Scalar QPSolver::get_function_value(GammaVector<Scalar> x){
	return this->get_function_value(x,false);
}

Scalar QPSolver::get_function_value(GammaVector<Scalar> x, bool use_gradient){
	Scalar fx = 0.0;

	if(use_gradient){
		/* use computed gradient in this->gs to compute function value */
		fx = 0.5*get_dot(this->g-this->b,x);
	} else {
		/* we have nothing - compute fx using full formula fx = 0.5*dot(A*x,x) - dot(b,x) */
		/* for safety - do not use any allocated vector */
		
		GammaVector<Scalar> Ax(this->get_T()*this->get_K());
		Scalar xAx, xb;

		get_Ax_laplace(Ax,x,this->eps_sqr);
		 
		xAx = get_dot(Ax,x);
		fx += 0.5*xAx;
		 
		xb = get_dot(x,this->b);
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
	
	oss << oss_spaces.str() << "-- QP SOLVER --";
	Message_info(oss.str());
	oss.str(""); oss.clear();

	oss << oss_spaces.str() << " - K = ";
	Message_info_value(oss.str(),this->get_K());
	oss.str(""); oss.clear();
	
	oss << oss_spaces.str() << " - T = ";
	Message_info_value(oss.str(),this->get_T());
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
		oss_values << this->x(k*T,(k+1)*T-1);
		Message_info_values(oss.str(),oss_values.str());	
		oss.str(""); oss.clear();
		oss_values.str(""); oss_values.clear();
	}

}

int QPSolver::get_T(){
	return this->T;
}

int QPSolver::get_K(){
	return this->K;
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
	return this->timer_projection.get_value_sum();
}


double QPSolver::get_time_matmult(){
	return this->timer_matmult.get_value_sum();
}

double QPSolver::get_time_dot(){
	return this->timer_dot.get_value_sum();
}

double QPSolver::get_time_update(){
	return this->timer_update.get_value_sum();
}

double QPSolver::get_time_total(){
	return this->timer_total.get_value_sum();
}

double QPSolver::get_time_stepsize(){
	return this->timer_stepsize.get_value_sum();
}

double QPSolver::get_time_fs(){
	return this->timer_fs.get_value_sum();
}

