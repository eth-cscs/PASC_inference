#include "external/petscvector/solver/entropysolvernewton.h"

#ifdef USE_DLIB

namespace pascinference {
namespace solver {


template<> EntropySolverNewton<PetscVector>::EntropySolverNewton(EntropyData<PetscVector> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;
	int K = entropydata->get_K();

	/* initial values */
	this->it_sums = new int [K];
	this->it_lasts = new int [K];
	this->itAxb_sums = new int [K];
	this->itAxb_lasts = new int [K];
	this->fxs = new double [K];
	this->gnorms = new double [K];

	set_value_array(K, this->it_sums, 0);
	set_value_array(K, this->it_lasts, 0);
	set_value_array(K, this->itAxb_sums, 0);
	set_value_array(K, this->itAxb_lasts, 0);
	set_value_array(K, this->fxs, std::numeric_limits<double>::max());
	set_value_array(K, this->gnorms, std::numeric_limits<double>::max());

	this->number_of_moments = EntropyData<PetscVector>::compute_number_of_moments(get_xdim(), get_Km());

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_compute_moments.restart();	
	this->timer_solve.restart();
	this->timer_Axb.restart();
	this->timer_update.restart();
	this->timer_g.restart();
	this->timer_H.restart();
	this->timer_fs.restart();
	this->timer_integrate.restart();

	/* prepare external content with PETSc-DLIB stuff */
	externalcontent = new ExternalContent();

	allocate_temp_vectors();

	LOG_FUNC_END	
}


template<> 
void EntropySolverNewton<PetscVector>::allocate_temp_vectors(){
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
	TRYCXX( VecSetSizes(integrals_Vec,entropydata->get_K()*(2*entropydata->get_Km()+1),PETSC_DECIDE) );
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
	this->delta = new GeneralVector<PetscVector>(*g);


	/* create Hessian matrix */
	TRYCXX( MatCreate(PETSC_COMM_SELF, &(externalcontent->H_petsc)) );
	TRYCXX( MatSetSizes(externalcontent->H_petsc,entropydata->get_Km(),entropydata->get_Km(),PETSC_DECIDE,PETSC_DECIDE) );
	#ifdef USE_CUDA
		TRYCXX( MatSetType(externalcontent->H_petsc,MATSEQAIJCUSPARSE) ); 
	#else
		TRYCXX( MatSetType(externalcontent->H_petsc,MATSEQAIJ) ); 
	#endif
	TRYCXX( MatSetFromOptions(externalcontent->H_petsc) );
	TRYCXX( MatMPIAIJSetPreallocation(externalcontent->H_petsc,entropydata->get_Km(),NULL,entropydata->get_Km()-1,NULL) ); 
	TRYCXX( MatSeqAIJSetPreallocation(externalcontent->H_petsc,entropydata->get_Km(),NULL) );

	TRYCXX( MatAssemblyBegin(externalcontent->H_petsc,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(externalcontent->H_petsc,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)(externalcontent->H_petsc),"Hessian matrix") );

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(x_power);
	free(x_power_gammak);
	free(moments_data);
	free(integrals);

	free(g);
	free(delta);

	TRYCXX( MatDestroy(&(externalcontent->H_petsc)) );

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start();
	 this->compute_moments_data();
	this->timer_compute_moments.stop();
	
//	coutMaster << "Moments: " << *moments << std::endl;

	this->timer_solve.start(); 
	int it; /* actual number of iterations */
	int itAxb, itAxb_one_newton_iteration; /* number of all KSP iterations, number of KSP iterations in one newton iteration */
	
	/* get dimensions */
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* get PETSc vecs - LOCAL */
	Vec moments_Vec = moments_data->get_vector();
	Vec x_Vec = entropydata->get_lambda()->get_vector(); /* x:=lambda is unknowns */
	Vec g_Vec = g->get_vector();
	Vec delta_Vec = delta->get_vector();
	Vec integrals_Vec = integrals->get_vector();

	/* stuff for clusters (subvectors) - LOCAL */
	IS k_is;
	IS integralsk_is;
	Vec xk_Vec;
	Vec momentsk_Vec;
	Vec integralsk_Vec;

	double *xk_arr;
	double *integralsk_arr;

	double gnorm, fx; /* norm(g),f(x) */
	double deltanorm = std::numeric_limits<double>::max(); /* norm(delta) - for debug */
	double deltanorm_old;

	/* postprocess */
	Vec g_inner_Vec;
	double gnorm_inner;

	/* through all clusters */
	for(int k = 0; k < K; k++){
		
		/* print iteration info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}
		
		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km, k*Km, 1, &k_is) ); /* Theta is LOCAL ! */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, Km+1, k*(Km+1), 1, &integralsk_is) ); 
	
		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(x_Vec, k_is, &xk_Vec) );
		TRYCXX( VecGetSubVector(moments_Vec, k_is, &momentsk_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* -------------- Newton algorithm (for k-th problem) ------------ */
		it = 0;
		itAxb = 0;
		
		/* compute integrals and gradient */
		this->timer_integrate.start();
		 TRYCXX( VecGetArray(xk_Vec,&xk_arr));
		 TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
 		  entropyintegration->compute(integralsk_arr, xk_arr, 2*Km+1);
		 TRYCXX( VecRestoreArray(xk_Vec,&xk_arr));
		 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
		this->timer_integrate.stop();
		this->timer_g.start();
		 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_g.stop();
		this->timer_fs.start();
		 fx = externalcontent->compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_fs.stop();

		while(it < this->maxit){
			/* compute stopping criteria - norm of gradient */
			TRYCXX( VecNorm(g_Vec, NORM_2, &gnorm) );
			if(gnorm < this->eps){
				break;
			}

			/* prepare Hessian matrix */
			this->timer_H.start();
			 externalcontent->compute_hessian(integralsk_Vec);
			this->timer_H.stop();
			
			/* ------ KSP Solver ----- */
			/* for solving H*delta=-g */
			TRYCXX( VecScale(g_Vec,-1.0) );
			
			/* Create linear solver context */
			TRYCXX( KSPCreate(PETSC_COMM_SELF,&(externalcontent->ksp)) ); // TODO: really every iteration?
			
			/* Set operators. Here the matrix that defines the linear system */
			/* also serves as the preconditioning matrix. */
			TRYCXX( KSPSetOperators(externalcontent->ksp,externalcontent->H_petsc,externalcontent->H_petsc) ); 
			
			/* precondition - maybe we can try LU factorization - then system will be solved in one iteration */
			TRYCXX( KSPGetPC(externalcontent->ksp,&(externalcontent->pc)) );
			TRYCXX( PCSetType(externalcontent->pc,PCJACOBI) );
			
			/* set stopping criteria */
			TRYCXX( KSPSetTolerances(externalcontent->ksp,this->eps_Axb,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT) ); // TODO: look onto all arguments

			/* some funny stuff can be loaded from petsc console parameters */
			TRYCXX( KSPSetFromOptions(externalcontent->ksp) );
			
			/* I think that these thing will use actual values in delta as initial approximation */
			TRYCXX( KSPSetInitialGuessNonzero(externalcontent->ksp,PETSC_TRUE) );
			
			/* Solve linear system */
			this->timer_Axb.start();
			 TRYCXX( KSPSolve(externalcontent->ksp,g_Vec,delta_Vec) );
			this->timer_Axb.stop();
			
			/* get some informations about solution process */
			TRYCXX( KSPGetIterationNumber(externalcontent->ksp,&itAxb_one_newton_iteration) );
			itAxb += itAxb_one_newton_iteration;
			
			/* print info about KSP */
			if(debug_print_Axb){
				TRYCXX( KSPView(externalcontent->ksp,PETSC_VIEWER_STDOUT_SELF) );
			}
			
			/* destroy KSP solver */
			TRYCXX( KSPDestroy(&(externalcontent->ksp)) ); 

			/* compute norm of inner error */
			if(debug_print_it){
				/* compute norm(H*delta-g) */
				TRYCXX( VecDuplicate(g_Vec,&g_inner_Vec) );
				TRYCXX( MatMult(externalcontent->H_petsc,delta_Vec,g_inner_Vec) );
				TRYCXX( VecAXPY(g_inner_Vec,-1.0,g_Vec) );
				
				TRYCXX( VecNorm(g_inner_Vec, NORM_2, &gnorm_inner) );
			}
			
			/* use solution from KSP to update newton iterations */
			this->timer_update.start();
			 TRYCXX( VecAXPY(xk_Vec,this->newton_coeff,delta_Vec)); /* x = x + delta; */
			this->timer_update.stop();
			
			/* recompute integrals, gradient, function value */
			this->timer_integrate.start();
			 TRYCXX( VecGetArray(xk_Vec,&xk_arr));
			 TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
			  entropyintegration->compute(integralsk_arr, xk_arr, 2*Km+1);
			 TRYCXX( VecRestoreArray(xk_Vec,&xk_arr));
			 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
			this->timer_integrate.stop();
			this->timer_g.start();
			 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_g.stop();
			this->timer_fs.start();
			 fx = externalcontent->compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_fs.stop();			
			
			it++;

			deltanorm_old = deltanorm;
			TRYCXX( VecNorm(delta_Vec, NORM_2, &deltanorm) );

			if(deltanorm > deltanorm_old){
				break; //TODO: hotfix
			}

			/* print progress of algorithm */
			if(debug_print_it){
				TRYCXX( VecNorm(delta_Vec, NORM_2, &deltanorm) );
				
				coutMaster << "\033[33m   it = \033[0m" << it;
				std::streamsize ss = std::cout.precision();
				coutMaster << ", \t\033[36mfx = \033[0m" << std::setprecision(17) << fx << std::setprecision(ss);
				coutMaster << ", \t\033[36mnorm(delta) = \033[0m" << std::setprecision(17) << deltanorm << std::setprecision(ss);
				coutMaster << ", \t\033[36mnorm(g_inner) = \033[0m" << gnorm_inner;
				coutMaster << ", \t\033[36mnorm(g_outer) = \033[0m" << gnorm << std::endl;

				/* log function value */
				LOG_FX(fx)
			}


			///* monitor - export values of stopping criteria */
			//if(this->monitor && GlobalManager.get_rank() == 0){
				////TODO: this could be done in a different way
				//std::ofstream myfile;
				//myfile.open("log/entropysolvernewton_monitor.m", std::fstream::in | std::fstream::out | std::fstream::app);

				//std::streamsize ss = myfile.precision();
				//myfile << std::setprecision(17);
			
				//myfile << "fx(" << it << ") = " << fx << "; "; 
				//myfile << "alpha_bb(" << it << ") = " << alpha_bb << "; ";
				//myfile << "gnorm(" << it << ") = " << gnorm << "; ";
				//myfile << std::endl;
			
				//myfile << std::setprecision(ss);
				//myfile.close();			
			//}

		}

		/* store number of iterations */
		this->it_lasts[k] = it;
		this->it_sums[k] += it;
		this->itAxb_lasts[k] = itAxb;
		this->itAxb_sums[k] += itAxb;
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
void EntropySolverNewton<PetscVector>::compute_moments_data() {
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
void EntropySolverNewton<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
	LOG_FUNC_BEGIN

	int T = entropydata->get_T();
	int Tlocal = entropydata->get_decomposition()->get_Tlocal();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	double F_;

	Vec lambda_Vec = entropydata->get_lambda()->get_vector();
	IS lambdak_is;
	Vec lambdak_Vec;
	double *lambdak_arr;

	Vec integrals_Vec = integrals->get_vector();
	IS integralsk_is;
	Vec integralsk_Vec;
	double *integralsk_arr;
	
	const double *x_arr;
	TRYCXX( VecGetArrayRead(entropydata->get_x()->get_vector(), &x_arr) );
	
	double *residuum_arr;
	TRYCXX( VecGetArray(residuum->get_vector(), &residuum_arr) );

	double mysum, x_power;
	for(int k=0;k<K;k++){

		/* get lambdak */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF,Km, k*Km, 1, &lambdak_is) );
		TRYCXX( VecGetSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( VecGetArray(lambdak_Vec, &lambdak_arr) );

		/* get integrals */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF,2*Km+1, k*(2*Km+1), 1, &integralsk_is) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );
		TRYCXX( VecGetArray(integralsk_Vec, &integralsk_arr) );

		/* compute int_X exp(-sum lambda_j x^j) dx for this cluster */
		/* = only integral with km=0 */
		entropyintegration->compute(integralsk_arr, lambdak_arr, 1);

		F_ = log(integralsk_arr[0]);

		for(int t=0;t<Tlocal;t++){
			/* compute sum_{j=1}^{Km} lambda_j*x^j */
			mysum = 0.0;
			x_power = x_arr[t]; /* x^1 */
			for(int km=0;km<Km;km++){
				mysum += lambdak_arr[km]*x_power;
				x_power *= x_arr[t]; /* x_power = x^(km+1) */
			}

			residuum_arr[t*K + k] = mysum + F_;
		}

		/* restore arrays */
		TRYCXX( VecRestoreArray(integralsk_Vec, &integralsk_arr) );
		TRYCXX( VecRestoreArray(lambdak_Vec, &lambdak_arr) );

		/* restore integralsk */
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );
		TRYCXX( ISDestroy(&integralsk_is) );

		/* restore lambdak */
		TRYCXX( VecRestoreSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( ISDestroy(&lambdak_is) );
	}

	/* restore arrays */
	TRYCXX( VecRestoreArray(residuum->get_vector(), &residuum_arr) );
	TRYCXX( VecRestoreArrayRead(entropydata->get_x()->get_vector(), &x_arr) );

		
	LOG_FUNC_END
}

/* ----------------------- external content */
void EntropySolverNewton<PetscVector>::ExternalContent::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN

	int Km;
	TRYCXX( VecGetSize(moments_Vec, &Km) );
    
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

void EntropySolverNewton<PetscVector>::ExternalContent::compute_hessian(Vec &integrals_Vec) {
	LOG_FUNC_BEGIN

	int Km;
	TRYCXX( MatGetSize(H_petsc, &Km, NULL) );
	    
    double *integrals_arr;

	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );

    /* fill Hessian matrix */
	for(int km1=0; km1 < Km; km1++){
		for(int km2=0; km2 < Km; km2++){
			double value =  integrals_arr[km1+km2+2]/integrals_arr[0] - integrals_arr[km1+1]*integrals_arr[km2+1]/(integrals_arr[0]*integrals_arr[0]);
			TRYCXX( MatSetValue(H_petsc, km1, km2, value, INSERT_VALUES) );
		}
	} 

	/* assemble matrix */
	TRYCXX( MatAssemblyBegin(H_petsc,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(H_petsc,MAT_FINAL_ASSEMBLY) );	

	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
	
	LOG_FUNC_END
}

double EntropySolverNewton<PetscVector>::ExternalContent::compute_function_value(Vec &lambda_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
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

template<> EntropySolverNewton<PetscVector>::ExternalContent * EntropySolverNewton<PetscVector>::get_externalcontent() const {
	return externalcontent;
}



}
} /* end namespace */

#endif
