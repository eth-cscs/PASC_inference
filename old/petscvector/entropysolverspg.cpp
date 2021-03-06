#include "external/petscvector/solver/entropysolverspg.h"

#ifdef USE_DLIB

namespace pascinference {
namespace solver {

template<>
void EntropySolverSPG<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	/* prepare auxiliary vectors */
	x_power = new GeneralVector<PetscVector>(*entropydata->get_x());
	x_power_gammak = new GeneralVector<PetscVector>(*entropydata->get_x());
	
	/* create aux vector for the computation of moments and integrals */
	Vec moments_Vec;
	Vec integrals_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&moments_Vec) );
	TRYCXX( VecCreate(PETSC_COMM_SELF,&integrals_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(moments_Vec, VECSEQCUDA));
		TRYCXX(VecSetType(integrals_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(moments_Vec, VECSEQ));
		TRYCXX(VecSetType(integrals_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*entropydata->get_Km(),PETSC_DECIDE) );
	TRYCXX( VecSetSizes(integrals_Vec,entropydata->get_K()*(entropydata->get_Km()+1),PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	TRYCXX( VecSetFromOptions(integrals_Vec) );
	this->moments_data = new GeneralVector<PetscVector>(moments_Vec);
	this->integrals = new GeneralVector<PetscVector>(integrals_Vec);

	/* SPG stuff - create vector to solve k-th problem (of size Km) */
	Vec g_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&g_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(g_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(g_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(g_Vec,entropydata->get_Km(),entropydata->get_Km()) );
	TRYCXX( VecSetFromOptions(g_Vec) );
	this->g = new GeneralVector<PetscVector>(g_Vec);
	this->y = new GeneralVector<PetscVector>(*g);
	this->s = new GeneralVector<PetscVector>(*g);

	LOG_FUNC_END
}

/* destroy temp_vectors */
template<>
void EntropySolverSPG<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(x_power);
	free(x_power_gammak);
	free(moments_data);
	free(integrals);

	free(g);
	free(y);
	free(s);

	LOG_FUNC_END
}

template<>
void EntropySolverSPG<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start();
	 this->compute_moments_data();
	this->timer_compute_moments.stop();
	
//	coutMaster << "Moments: " << *moments << std::endl;

	this->timer_solve.start(); 
	int it, itgll; /* actual number of iterations */
	SPG_fs fs(this->m); /* store function values for generalized Armijo condition (GLL) */

	/* get dimensions */
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* get PETSc vecs - LOCAL */
	Vec moments_Vec = moments_data->get_vector();
	Vec x_Vec = entropydata->get_lambda()->get_vector(); /* x:=lambda is unknowns */
	Vec g_Vec = g->get_vector();
	Vec y_Vec = y->get_vector();
	Vec s_Vec = s->get_vector();
	Vec integrals_Vec = integrals->get_vector();

	/* stuff for clusters (subvectors) - LOCAL */
	IS k_is;
	IS integralsk_is;
	Vec xk_Vec;
	Vec momentsk_Vec;
	Vec integralsk_Vec;

	double alpha_bb, gnorm, fx, fx_old, fx_max, delta, beta, beta_temp; /* BB & GLL stuff */
	int itgll_temp;
	double sTs, sTy;

	/* through all clusters */
	for(int k = 0; k < K; k++){
		
		coutMaster << "k=" << k << std::endl;
		
		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km, k*Km, 1, &k_is) ); /* Theta is LOCAL ! */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km+1, k*(Km+1), 1, &integralsk_is) ); 
	
		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(x_Vec, k_is, &xk_Vec) );
		TRYCXX( VecGetSubVector(moments_Vec, k_is, &momentsk_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* -------------- SPG algorithm (for k-th problem) ------------ */
		it = 0;
		itgll = 0;
		
		/* compute integrals and gradient */
		this->timer_integrate.start();
		 compute_integrals(integralsk_Vec, xk_Vec, true);
		this->timer_integrate.stop();
		this->timer_g.start();
		 compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_g.stop();
		this->timer_fs.start();
		 fx = compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
		 fs.init(fx);
		this->timer_fs.stop();

		/* set initial step-size */
		alpha_bb = this->alphainit;

		while(it < this->maxit){
			/* compute stopping criteria - norm of gradient */
			TRYCXX( VecNorm(g_Vec, NORM_2, &gnorm) );
			if(gnorm < this->eps){
				break;
			}

			/* store previous iteration */
			fx_old = fx;
			TRYCXX( VecCopy(xk_Vec,s_Vec) ); /* s := x_old; */
			TRYCXX( VecCopy(g_Vec,y_Vec) ); /* y := g_old; */
		
			/* update */
			this->timer_update.start();
			 TRYCXX( VecAXPY(xk_Vec, -alpha_bb, y_Vec) ); /* x = x - alpha_bb*g */
			this->timer_update.stop();

			/* we have new approximation => compute integrals and gradient */
			this->timer_integrate.start();
			 compute_integrals(integralsk_Vec, xk_Vec, false);
			this->timer_integrate.stop();
			this->timer_fs.start();
			 fx = compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_fs.stop();

			/* ------------- modify approximation using GLL ------------- */
			/* fx_old := f(x_old)  = f(x0); */
			/* s      := x_old     = x0; */
			/* y      := g_old     = g0; */ 
			/* fx     := f(x_temp); */
			/* x      := x_temp; */
			/* g      := g_temp; */
			
				/* delta = dot(g0,d) */
				this->timer_dot.start();
				 TRYCXX( VecDot(y_Vec,y_Vec,&delta) );
				this->timer_dot.stop();
				delta = -alpha_bb*delta;

				/* fx_max = max(fs) */
				this->timer_fs.start();
				 fx_max = fs.get_max();
				this->timer_fs.stop();

				/* initial cutting parameter */
				beta = 1.0;
				itgll_temp = 0;
				
				/* while f_temp > f_max + gamma*beta*delta */
				while(fx > fx_max + this->gamma*beta*delta && itgll_temp < this->maxit_gll ){
					/* beta_temp = -0.5*beta^2*delta/(f_temp - f_0 - beta*delta); */
					beta_temp = -0.5*beta*beta*delta/(fx - fx_old - beta*delta);
					
					if(beta_temp >= this->sigma1 && beta_temp <= this->sigma2*beta){
						beta = beta_temp;
					} else {
						beta = 0.9*beta;
					}
					
					/* update approximation */
					this->timer_update.start();
					 TRYCXX( VecCopy(s_Vec, xk_Vec) ); /* x = x_old */
					 TRYCXX( VecAXPY(xk_Vec, (-1)*alpha_bb*beta, y_Vec) ); /* x = x + beta*g */
					this->timer_update.stop();

					/* we have new approximation => compute integrals */
					this->timer_integrate.start();
					 compute_integrals(integralsk_Vec, xk_Vec, false);
					this->timer_integrate.stop();
					this->timer_fs.start();
					 fx = compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
					this->timer_fs.stop();
					
					itgll_temp++;
				}
				
				/* update gll iterator */
				itgll += itgll_temp;

			/* ------------- end of GLL ------------- */

			/* store last computed function value */
			this->timer_fs.start();
			 fs.update(fx);
			this->timer_fs.stop();

			/* compute all integrals and gradient in new approximation */
			this->timer_integrate.start();
			 compute_integrals(integralsk_Vec, xk_Vec, true);
			this->timer_integrate.stop();
			this->timer_g.start();
			 compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_g.stop();


			/* compute s = x - x_old; y = g - g_old; */
			this->timer_update.start();
			 TRYCXX( VecAYPX(s_Vec, -1.0, xk_Vec) );
			 TRYCXX( VecAYPX(y_Vec, -1.0, g_Vec) );
			this->timer_update.stop();

			/* compute new bb-step */
			this->timer_dot.start();
			 TRYCXX( VecDot(s_Vec, s_Vec, &sTs) );
			 TRYCXX( VecDot(s_Vec, y_Vec, &sTy) );
			this->timer_dot.stop();

			alpha_bb = sTs/sTy;
//			alpha_bb = sTy/sTs;

			it++;

			/* print progress of algorithm */
			if(debug_print_it){
				coutMaster << "\033[33m   it = \033[0m" << it;
			
				std::streamsize ss = std::cout.precision();
				coutMaster << ", \t\033[36mfx = \033[0m" << std::setprecision(17) << fx << std::setprecision(ss);
				coutMaster << ", \t\033[36malpha_bb = \033[0m" << std::setprecision(17) << alpha_bb << std::setprecision(ss);
				coutMaster << ", \t\033[36mbeta_gll = \033[0m" << std::setprecision(17) << beta << std::setprecision(ss);
				coutMaster << ", \t\033[36mgnorm = \033[0m" << gnorm << std::endl;

				/* log function value */
				LOG_FX(fx)
			}

			/* monitor - export values of stopping criteria */
			if(this->monitor && GlobalManager.get_rank() == 0){
				//TODO: this could be done in a different way
				std::ofstream myfile;
				myfile.open("log/entropysolverspg_monitor.m", std::fstream::in | std::fstream::out | std::fstream::app);

				std::streamsize ss = myfile.precision();
				myfile << std::setprecision(17);
			
				myfile << "fx(" << it << ") = " << fx << "; "; 
				myfile << "alpha_bb(" << it << ") = " << alpha_bb << "; ";
				myfile << "gnorm(" << it << ") = " << gnorm << "; ";
				myfile << std::endl;
			
				myfile << std::setprecision(ss);
				myfile.close();			
			}

		}

		/* store number of iterations */
		this->it_lasts[k] = it;
		this->it_sums[k] += it;
		this->itgll_lasts[k] = itgll;
		this->itgll_sums[k] += itgll;
		this->fxs[k] = fx;
		this->gnorms[k] = gnorm;

		/* -------------- end of SPG algorithm ------------ */

		/* restore subvectors */
		TRYCXX( VecRestoreSubVector(x_Vec, k_is, &xk_Vec) );
		TRYCXX( VecRestoreSubVector(moments_Vec, k_is, &momentsk_Vec) );
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* destroy index set */
		TRYCXX( ISDestroy(&k_is) );
		TRYCXX( ISDestroy(&integralsk_is) );	

	} /* endfor through clusters */


	this->timer_solve.stop(); 

	LOG_FUNC_END
}

template<>
void EntropySolverSPG<PetscVector>::compute_moments_data() {
	LOG_FUNC_BEGIN

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec x_power_Vec = x_power->get_vector();
	Vec x_power_gammak_Vec = x_power_gammak->get_vector();
	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments_data->get_vector();
	
	Vec gammak_Vec;
	IS gammak_is;
	
	TRYCXX( VecCopy(x_Vec,x_power_Vec) ); /* x^1 */
	
	double *moments_arr, mysum, gammaksum;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
	for(int km=0; km < entropydata->get_Km(); km++){
		
		for(int k=0;k<entropydata->get_K();k++){
			/* get gammak */
			this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

			/* compute x_power_gammak */
			TRYCXX( VecPointwiseMult(x_power_gammak_Vec, gammak_Vec, x_power_Vec) ); /* x_power_gammak = x_power.*gammak */

			/* compute gammaksum */
			TRYCXX( VecSum(gammak_Vec, &gammaksum) );
			TRYCXX( VecSum(x_power_gammak_Vec, &mysum) );

			/* store computed moment */
			if(gammaksum != 0){
				moments_arr[k*this->entropydata->get_Km() + km] = mysum/gammaksum;
			} else {
				moments_arr[k*this->entropydata->get_Km() + km] = 0.0;
			}
	
			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );	
		}
		
		TRYCXX( VecPointwiseMult(x_power_Vec, x_Vec, x_power_Vec) ); /* x_power = x_power.*x */
	}
	
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	LOG_FUNC_END
}

template<>
void EntropySolverSPG<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
	LOG_FUNC_BEGIN

	int T = entropydata->get_T();
	int Tlocal = entropydata->get_decomposition()->get_Tlocal();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* lambda vector for Dlib integration */
	column_vector lambda(Km);
    auto mom_function = [&](double x)->double { return gg(x, 0, lambda);};
    double F_;

	/* update gamma_solver data - prepare new linear term */
	/* theta includes all moments */
	const double *lambda_arr;
	TRYCXX( VecGetArrayRead(entropydata->get_lambda()->get_vector(), &lambda_arr) );
	
	const double *x_arr;
	TRYCXX( VecGetArrayRead(entropydata->get_x()->get_vector(), &x_arr) );
	
	double *residuum_arr;
	TRYCXX( VecGetArray(residuum->get_vector(), &residuum_arr) );

	double mysum, x_power;
	for(int k=0;k<K;k++){
		/* compute int_X exp(-sum lambda_j x^j) dx for this cluster
		/* from arr to Dlib-vec */
		for(int km=0;km<Km;km++){
			lambda(km) = lambda_arr[k*Km+km];
		}
		F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
		F_ = log(F_);

		for(int t=0;t<Tlocal;t++){
			/* compute sum_{j=1}^{Km} lambda_j*x^j */
			mysum = 0.0;
			x_power = x_arr[t]; /* x^1 */
			for(int km=0;km<Km;km++){
				mysum += lambda_arr[k*Km+km]*x_power;
				x_power *= x_arr[t]; /* x_power = x^(km+1) */
			}

			residuum_arr[t*K + k] = mysum + F_;
		}
	}

	/* coeffs of A_shared are updated via computation of Theta :) */

	/* restore arrays */
	TRYCXX( VecRestoreArray(residuum->get_vector(), &residuum_arr) );
	TRYCXX( VecRestoreArrayRead(entropydata->get_x()->get_vector(), &x_arr) );
	TRYCXX( VecRestoreArrayRead(entropydata->get_lambda()->get_vector(), &lambda_arr) );
		
	LOG_FUNC_END
}

template<>
void EntropySolverSPG<PetscVector>::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN

	int Km = entropydata->get_Km();
    
    double *g_arr;
    double *integrals_arr;
    double *moments_arr;

	TRYCXX( VecGetArray(g_Vec, &g_arr) );
	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
    
    /* compute gradient */
    for (int km = 0; km < Km; km++){
        g_arr[km] = moments_arr[km] - integrals_arr[km+1]/integrals_arr[0];
	} 

	TRYCXX( VecRestoreArray(g_Vec, &g_arr) );
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );
	
	LOG_FUNC_END
}

template<>
double EntropySolverSPG<PetscVector>::compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN
	
	double f, momTlambda;
    double *integrals_arr;

	TRYCXX( VecDot(moments_Vec, lambda_Vec, &momTlambda) );

	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );
	f = log(integrals_arr[0]) + momTlambda;
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
		
	LOG_FUNC_END

	return f;
}

template<>
void EntropySolverSPG<PetscVector>::compute_integrals(Vec &integrals_Vec, Vec &lambda_Vec, bool compute_all) {
	LOG_FUNC_BEGIN

	int Km = entropydata->get_Km();
	int Km_int; /* number of computed integrals */
	if(compute_all){
		Km_int = Km;
	} else {
		Km_int = 0;
	}

	double *integrals_arr;
	double *lambda_arr;
	TRYCXX( VecGetArray(lambda_Vec, &lambda_arr) );
	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );

	//// TODO: temp
    column_vector lambda_Dlib(Km);
    for(int km=0;km<Km;km++){
		lambda_Dlib(km) = lambda_arr[km];
	}

    /* compute integrals */
    for(int km = 0; km<=Km_int;km++){
		auto mom_function = [&](double x)->double { return gg(x, km, lambda_Dlib);};
		integrals_arr[km] = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
	}

	TRYCXX( VecRestoreArray(lambda_Vec, &lambda_arr) );
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );

	LOG_FUNC_END
}



}
} /* end namespace */

#endif
