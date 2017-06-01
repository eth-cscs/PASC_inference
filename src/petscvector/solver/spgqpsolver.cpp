#include "external/petscvector/solver/spgqpsolver.h"

namespace pascinference {
namespace solver {

template<>
std::string SPGQPSolver<PetscVector>::get_name() const {
	return "SPGQPSolver for PETSc"; /* better to see than simple "SPGQPSolver<PetscVector>" */
}

/* prepare temp_vectors */
template<>
void SPGQPSolver<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<PetscVector> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<PetscVector>(*pattern);
	d = new GeneralVector<PetscVector>(*pattern);
	Ad = new GeneralVector<PetscVector>(*pattern);	
	temp = new GeneralVector<PetscVector>(*pattern);	

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();

	/* for Mdot */
	TRYCXX( PetscMalloc1(3,&(Mdots_val)) );
	TRYCXX( PetscMalloc1(3,&(externalcontent->Mdots_vec)) );

	externalcontent->Mdots_vec[0] = d->get_vector();
	externalcontent->Mdots_vec[1] = Ad->get_vector();
	externalcontent->Mdots_vec[2] = g->get_vector();

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<>
void SPGQPSolver<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(g);
	free(d);
	free(Ad);
	free(temp);
	
	TRYCXX( PetscFree(Mdots_val) );
	TRYCXX( PetscFree(externalcontent->Mdots_vec) );
	
	LOG_FUNC_END
}

/* solve the problem */
template<>
void SPGQPSolver<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	/* get Petsc objects from general */
	BlockGraphSparseMatrix<PetscVector> *Abgs = dynamic_cast<BlockGraphSparseMatrix<PetscVector> *>(qpdata->get_A());
	GeneralVector<PetscVector> *b_p = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_b());
	GeneralVector<PetscVector> *x_p = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x());
	GeneralVector<PetscVector> *x0_p = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x0());
	GeneralVector<PetscVector> *g_p = dynamic_cast<GeneralVector<PetscVector> *>(this->g);
	GeneralVector<PetscVector> *d_p = dynamic_cast<GeneralVector<PetscVector> *>(this->d);
	GeneralVector<PetscVector> *Ad_p = dynamic_cast<GeneralVector<PetscVector> *>(this->Ad);

	Mat A_Mat = Abgs->get_externalcontent()->A_petsc; 
	Vec b_Vec = b_p->get_vector();
	Vec x_Vec = x_p->get_vector();
	Vec x0_Vec = x0_p->get_vector();
	Vec g_Vec = g_p->get_vector();
	Vec d_Vec = d_p->get_vector();
	Vec Ad_Vec = Ad_p->get_vector();

	#ifdef USE_CUDA
		/* make sure that we are computing on GPU */
		TRYCXX( VecCUDACopyToGPU(b_Vec) );	
		TRYCXX( VecCUDACopyToGPU(x_Vec) );	
		TRYCXX( VecCUDACopyToGPU(x0_Vec) );	
		TRYCXX( VecCUDACopyToGPU(g_Vec) );	
		TRYCXX( VecCUDACopyToGPU(d_Vec) );	
		TRYCXX( VecCUDACopyToGPU(Ad_Vec) );	
	#endif
	
	allbarrier<PetscVector>();

	this->timer_solve.start(); /* stop this timer in the end of solution */

	int it = 0; /* number of iterations */
	int hessmult = 0; /* number of hessian multiplications */

	double fx; /* function value */
	double fx_old; /* f(x_{it - 1}) */
	SPG_fs fs(this->m); /* store function values for generalized Armijo condition */
	double fx_max; /* max(fs) */
	double xi, beta_bar, beta_hat, beta; /* for Armijo condition */
	double dd; /* dot(d,d) */
	double gd; /* dot(g,d) */
	double dAd; /* dot(Ad,d) */
	double alpha_bb; /* BB step-size */
	double normb; /* norm of linear term used in stopping criteria */

	/* compute normb for stopping criteria */
	TRYCXX( VecNorm(b_Vec, NORM_2, &normb) );
	allbarrier<PetscVector>();

	/* initial step-size */
	alpha_bb = this->alphainit;

	//TODO: temp!
	/* x = x0; set approximation as initial */
	TRYCXX( VecCopy(x0_Vec, x0_Vec) );
	allbarrier<PetscVector>();

	this->timer_projection.start();
	 qpdata->get_feasibleset()->project(*x_p); /* project initial approximation to feasible set */
 	 allbarrier<PetscVector>();
	this->timer_projection.stop();

	/* compute gradient, g = A*x-b */
	this->timer_matmult.start();
	 TRYCXX( MatMult(A_Mat, x_Vec, g_Vec) );
	 TRYCXX( VecScale(g_Vec, Abgs->get_coeff()) );
	 allbarrier<PetscVector>();
	 hessmult += 1; /* there was muliplication by A */
	this->timer_matmult.stop();

	TRYCXX( VecAXPY(g_Vec, -1.0, b_Vec) );
	allbarrier<PetscVector>();

	/* initialize fs */
	this->timer_fs.start();
	 fx = get_fx();
	 fx_old = std::numeric_limits<double>::max();
	 this->fx = fx;
	 fs.init(fx);
	this->timer_fs.stop();

	/* main cycle */
	while(it < this->maxit){
		/* increase iteration counter */
		it += 1;

		/* d = x - alpha_bb*g, see next step, it will be d = P(x - alpha_bb*g) - x */
		this->timer_update.start();
		 TRYCXX( VecCopy(x_Vec, d_Vec));
		 TRYCXX( VecAXPY(d_Vec, -alpha_bb, g_Vec) );
		 allbarrier<PetscVector>();
		this->timer_update.stop();

		/* d = P(d) */
		this->timer_projection.start();
		 qpdata->get_feasibleset()->project(*d_p);
		this->timer_projection.stop();

		/* d = d - x */
		this->timer_update.start();
		 TRYCXX( VecAXPY(d_Vec, -1.0, x_Vec) );
		 allbarrier<PetscVector>();
		this->timer_update.stop();

		/* Ad = A*d */
		this->timer_matmult.start();
		 TRYCXX( MatMult(A_Mat, d_Vec, Ad_Vec) );
		 TRYCXX( VecScale(Ad_Vec, Abgs->get_coeff()) );
		 allbarrier<PetscVector>();
		 hessmult += 1;
		this->timer_matmult.stop();

		this->timer_dot.start();
		 compute_dots(&dd, &dAd, &gd);
		this->timer_dot.stop();

		/* fx_max = max(fs) */
		this->timer_fs.start();
		 fx_max = fs.get_max();
		 allbarrier<PetscVector>();
		this->timer_fs.stop();
		
		/* compute step-size from A-condition */
		this->timer_stepsize.start();
		 xi = (fx_max - fx)/dAd;
		 beta_bar = -gd/dAd;
		 beta_hat = this->gamma*beta_bar + PetscSqrtReal(this->gamma*this->gamma*beta_bar*beta_bar + 2*xi);

		 /* beta = max(sigma1,min(sigma2,beta_hat)) */
		 if(beta_hat < this->sigma1){
			 beta_hat = this->sigma1;
		 }
		 
		 if(beta_hat < this->sigma2){
			beta = beta_hat;
		 } else {
			beta = this->sigma2;
		 }
		this->timer_stepsize.stop();

		/* x = x + beta*d; g = g + beta*Ad; update approximation and gradient */
		this->timer_update.start();
		 TRYCXX( VecAXPY(x_Vec, beta, d_Vec) );
		 TRYCXX( VecAXPY(g_Vec, beta, Ad_Vec) );
		 allbarrier<PetscVector>();
		this->timer_update.stop();

		/* compute new function value using gradient and update fs list */
		this->timer_fs.start();
		 fx_old = fx;
//		 fx = get_fx(fx_old,beta,gd,dAd);
		 fx = get_fx();
		 fs.update(fx);
		 allbarrier<PetscVector>();
		this->timer_fs.stop();

		/* update BB step-size */
		this->timer_stepsize.start();
		 alpha_bb = dd/dAd;
		this->timer_stepsize.stop();

		this->gP = dd;

		/* print progress of algorithm */
		if(debug_print_it){
			coutMaster << "\033[33m   it = \033[0m" << it;
			
			std::streamsize ss = std::cout.precision();
			coutMaster << ", \t\033[36mfx = \033[0m" << std::setprecision(17) << fx << std::setprecision(ss);

			coutMaster << ", \t\033[36mgP = \033[0m" << this->gP;
			coutMaster << ", \t\033[36mdd = \033[0m" << dd << std::endl;

			/* log function value */
			LOG_FX(fx)

		}
		
		/* print qpdata */
		if(debug_print_vectors){
//			coutMaster << "x: " << x << std::endl;
//			coutMaster << "d: " << d << std::endl;
//			coutMaster << "g: " << g << std::endl;
//			coutMaster << "Ad: " << Ad << std::endl;
		}

		if(debug_print_scalars){
			coutMaster << "\033[36m    alpha_bb = \033[0m" << alpha_bb << ",";
			coutMaster << "\033[36m dAd = \033[0m" << dAd << ",";
			coutMaster << "\033[36m gd = \033[0m" << gd << std::endl;
			
			coutMaster << "\033[36m    fx = \033[0m" << fx << ",";
			coutMaster << "\033[36m fx_max = \033[0m" << fx_max << ",";
			coutMaster << "\033[36m xi = \033[0m" << xi << std::endl;
			
			coutMaster << "\033[36m    beta_bar = \033[0m" << beta_bar << ",";
			coutMaster << "\033[36m beta_hat = \033[0m" << beta_hat << ",";
			coutMaster << "\033[36m beta = \033[0m" << beta << std::endl;
			
		}

		/* stopping criteria */
		if( this->stop_difff && abs(fx - fx_old) < this->eps){
			break;
		}
		if(this->stop_normgp && dd < this->eps){
			break;
		}
		if(this->stop_normgp_normb && dd < this->eps*normb){
			break;
		}
		if(this->stop_Anormgp && dAd < this->eps){
			break;
		}
		if(this->stop_Anormgp_normb && dAd < this->eps*normb){
			break;
		}
		
		/* monitor - export values of stopping criteria */
		if(this->monitor && GlobalManager.get_rank() == 0){
			//TODO: this could be done in a different way
			std::ofstream myfile;
			myfile.open("log/spgqpsolver_monitor.m", std::fstream::in | std::fstream::out | std::fstream::app);

			std::streamsize ss = myfile.precision();
			myfile << std::setprecision(17);
			
			myfile << "fx(" << it << ") = " << fx << "; "; 
			myfile << "alpha_bb(" << it << ") = " << alpha_bb << "; ";
			myfile << "beta(" << it << ") = " << beta << "; ";
			myfile << "norm_difff(" << it << ") = " << abs(fx - fx_old) << "; ";
			myfile << "norm_gp(" << it << ") = " << dd << "; ";
			myfile << "norm_Agp(" << it << ") = " << dAd << "; ";
			myfile << std::endl;
			
			myfile << std::setprecision(ss);
			myfile.close();			
		}
		
	} /* main cycle end */

	this->it_sum += it;
	this->hessmult_sum += hessmult;
	this->it_last = it;
	this->hessmult_last = hessmult;

	this->fx = fx;
	this->timer_solve.stop();

	/* write info to log file */
	LOG_IT(it)
	LOG_FX(fx)

	LOG_FUNC_END
}

/* compute function value using inner *x and already computed *g */
template<>
double SPGQPSolver<PetscVector>::get_fx() const {
	LOG_FUNC_BEGIN
	
	double fx = std::numeric_limits<double>::max();

	/* use computed gradient in this->g to compute function value */
	double tempt; 

	/* get PETSc specific stuff from general */
	GeneralVector<PetscVector> *b_p = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_b());
	GeneralVector<PetscVector> *x_p = dynamic_cast<GeneralVector<PetscVector> *>(qpdata->get_x());
	GeneralVector<PetscVector> *g_p = dynamic_cast<GeneralVector<PetscVector> *>(this->g);
	GeneralVector<PetscVector> *temp_p = dynamic_cast<GeneralVector<PetscVector> *>(this->temp);

	Vec b_Vec = b_p->get_vector();
	Vec x_Vec = x_p->get_vector();
	Vec g_Vec = g_p->get_vector();
	Vec temp_Vec = temp_p->get_vector();

	/* temp = g - b; tempt = dot(temp,x); */
	TRYCXX( VecCopy(g_Vec,temp_Vec) );
	TRYCXX( VecAXPY(temp_Vec, -1.0, b_Vec));
	TRYCXX( VecDot(x_Vec,temp_Vec,&tempt));
	allbarrier<PetscVector>();

	fx = 0.5*tempt;

	LOG_FUNC_END
	return fx;	
}

template<>
void SPGQPSolver<PetscVector>::compute_dots(double *dd, double *dAd, double *gd) const {
	LOG_FUNC_BEGIN

	//TODO: PETSC Mdot is not working in GPU, there is a hotfix:
//	TRYCXX( VecMDot( Mdots_vec[0], 3, Mdots_vec, Mdots_val) );
//	*dd = Mdots_val[0];
//	*dAd = Mdots_val[1];
//	*gd = Mdots_val[2];

	TRYCXX( VecDot( externalcontent->Mdots_vec[0], externalcontent->Mdots_vec[0], dd) );
	allbarrier<PetscVector>();

	TRYCXX( VecDot( externalcontent->Mdots_vec[0], externalcontent->Mdots_vec[1], dAd) );
	allbarrier<PetscVector>();

	TRYCXX( VecDot( externalcontent->Mdots_vec[0], externalcontent->Mdots_vec[2], gd) );
	allbarrier<PetscVector>();

	LOG_FUNC_END
}

template<> 
SPGQPSolver<PetscVector>::ExternalContent * SPGQPSolver<PetscVector>::get_externalcontent() const {
	return this->externalcontent;	
}


}
}
