#include "external/petscvector/solver/entropysolvernewton.h"

//TODO: uncomment, or maybe remove?
//#ifdef USE_DLIB

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

	/* set initial Theta equal to 1 */
	Vec lambda_Vec = entropydata->get_lambda()->get_vector();
	TRYCXX( VecSet(lambda_Vec, 1.0) );

	LOG_FUNC_END
}

template<>
EntropySolverNewton<PetscVector>::~EntropySolverNewton(){
	LOG_FUNC_BEGIN

	free(this->it_sums);
	free(this->it_lasts);
	free(this->itAxb_sums);
	free(this->itAxb_lasts);
	free(this->fxs);
	free(this->gnorms);

	/* free temp vectors */
	free_temp_vectors();

	/* free tool for integration */
	if(this->entropyintegration){
		free(this->entropyintegration);
	}

	/* destroy external content */
	free(externalcontent);

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	int number_of_moments = get_number_of_moments();
	int n = number_of_moments-1;
	int number_of_integrals = 1 + n + (int)(0.5*n*(n+1));

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
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*number_of_moments,PETSC_DECIDE) );
	TRYCXX( VecSetSizes(integrals_Vec,entropydata->get_K()*number_of_integrals,PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	TRYCXX( VecSetFromOptions(integrals_Vec) );
	this->moments = new GeneralVector<PetscVector>(moments_Vec);
	this->integrals = new GeneralVector<PetscVector>(integrals_Vec);

	/* SPG stuff - create vector to solve k-th problem (of size Km) */
	Vec g_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&g_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(g_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(g_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(g_Vec,n,n) );
	TRYCXX( VecSetFromOptions(g_Vec) );
	this->g = new GeneralVector<PetscVector>(g_Vec);
	this->delta = new GeneralVector<PetscVector>(*g);


	/* create Hessian matrix */
	TRYCXX( MatCreate(PETSC_COMM_SELF, &(externalcontent->H_petsc)) );
	TRYCXX( MatSetSizes(externalcontent->H_petsc,n,n,PETSC_DECIDE,PETSC_DECIDE) );
	#ifdef USE_CUDA
		TRYCXX( MatSetType(externalcontent->H_petsc,MATSEQAIJCUSPARSE) );
	#else
		TRYCXX( MatSetType(externalcontent->H_petsc,MATSEQAIJ) );
	#endif
	TRYCXX( MatSetFromOptions(externalcontent->H_petsc) );
	TRYCXX( MatMPIAIJSetPreallocation(externalcontent->H_petsc,n,NULL,n-1,NULL) );
	TRYCXX( MatSeqAIJSetPreallocation(externalcontent->H_petsc,n,NULL) );

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
	free(moments);
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
	 this->compute_moments();
	this->timer_compute_moments.stop();

//	coutMaster << "Moments: " << *moments << std::endl;

	this->timer_solve.start();
	int it; /* actual number of iterations */
	int itAxb, itAxb_one_newton_iteration; /* number of all KSP iterations, number of KSP iterations in one newton iteration */

	/* get dimensions */
	int number_of_moments = entropydata->get_number_of_moments();
	int n = number_of_moments-1;
	int number_of_integrals = entropyintegration->get_number_of_integrals();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* get PETSc vecs - LOCAL */
	Vec moments_Vec = moments->get_vector();
	Vec lambda_Vec = entropydata->get_lambda()->get_vector(); /* x:=lambda is unknowns */
	Vec g_Vec = g->get_vector();
	Vec delta_Vec = delta->get_vector();
	Vec integrals_Vec = integrals->get_vector();

    /* restart solver - initial approximation */
	TRYCXX( VecSet(lambda_Vec, 1.0) );

	/* stuff for clusters (subvectors) - LOCAL */
	IS lambdak_is;
	IS momentsk_is;
	IS integralsk_is;
	Vec lambdak_Vec;
	Vec momentsk_Vec;
	Vec integralsk_Vec;

	double *lambdak_arr;
	double *integralsk_arr;

	double gnorm, fx; /* norm(g),f(x) */
	double deltanorm = std::numeric_limits<double>::max(); /* norm(delta) - for debug */

	/* postprocess */
	Vec g_inner_Vec;
	double gnorm_inner = std::numeric_limits<double>::max();

	/* through all clusters */
	for(int k = 0; k < K; k++){

		/* print iteration info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}

		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, n, k*number_of_moments+1, 1, &lambdak_is) );
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, n, k*number_of_moments+1, 1, &momentsk_is) );
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, number_of_integrals, k*number_of_integrals, 1, &integralsk_is) );

		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( VecGetSubVector(moments_Vec, momentsk_is, &momentsk_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* -------------- Newton algorithm (for k-th problem) ------------ */
		it = 0;
		itAxb = 0;

		/* compute integrals and gradient */
		this->timer_integrate.start();
		 TRYCXX( VecGetArray(lambdak_Vec,&lambdak_arr));
		 TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
 		  entropyintegration->compute(integralsk_arr, lambdak_arr, number_of_integrals);
		 TRYCXX( VecRestoreArray(lambdak_Vec,&lambdak_arr));
		 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
		this->timer_integrate.stop();
		this->timer_g.start();
		 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_g.stop();
		this->timer_fs.start();
		 fx = externalcontent->compute_function_value(lambdak_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_fs.stop();

		while(it < this->maxit){
			/* compute stopping criteria - norm of gradient */
			TRYCXX( VecNorm(g_Vec, NORM_2, &gnorm) );

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

            /* use stopping criteria */
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
//			TRYCXX( PCSetType(externalcontent->pc,PCJACOBI) );
			TRYCXX( PCSetType(externalcontent->pc,PCNONE) );

			/* set stopping criteria */
			TRYCXX( KSPSetTolerances(externalcontent->ksp,this->eps_Axb,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT) ); // TODO: look onto all arguments

			/* some funny stuff can be loaded from petsc console parameters */
			TRYCXX( KSPSetFromOptions(externalcontent->ksp) );

			/* I think that these thing will use actual values in delta as initial approximation */
			TRYCXX( KSPSetInitialGuessNonzero(externalcontent->ksp,PETSC_TRUE) );

            /* before solving the system, print content */
            if(debug_print_vectors){
            	coutMaster << "lambdak : " << std::endl;
                TRYCXX( VecView(lambdak_Vec, PETSC_VIEWER_STDOUT_WORLD) );
                coutMaster << "moments : " << std::endl;
                TRYCXX( VecView(momentsk_Vec, PETSC_VIEWER_STDOUT_WORLD) );
                coutMaster << "integrals : " << std::endl;
                TRYCXX( VecView(integralsk_Vec, PETSC_VIEWER_STDOUT_WORLD) );
            	coutMaster << "gradient : " << std::endl;
                TRYCXX( VecView(g_Vec, PETSC_VIEWER_STDOUT_WORLD) );
            	coutMaster << "Hessian : " << std::endl;
                TRYCXX( MatView(externalcontent->H_petsc, PETSC_VIEWER_STDOUT_WORLD) );
            }


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
			 TRYCXX( VecAXPY(lambdak_Vec,this->newton_coeff,delta_Vec)); /* x = x + delta; */
			this->timer_update.stop();

			/* recompute integrals, gradient, function value */
			this->timer_integrate.start();
			 TRYCXX( VecGetArray(lambdak_Vec,&lambdak_arr));
			 TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
			  entropyintegration->compute(integralsk_arr, lambdak_arr, number_of_integrals);
			 TRYCXX( VecRestoreArray(lambdak_Vec,&lambdak_arr));
			 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
			this->timer_integrate.stop();
			this->timer_g.start();
			 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_g.stop();
			this->timer_fs.start();
			 fx = externalcontent->compute_function_value(lambdak_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_fs.stop();

			it++;

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

		/* -------------- end of Newton algorithm ------------ */

		/* restore subvectors */
		TRYCXX( VecRestoreSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( VecRestoreSubVector(moments_Vec, momentsk_is, &momentsk_Vec) );
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* destroy index set */
		TRYCXX( ISDestroy(&lambdak_is) );
		TRYCXX( ISDestroy(&momentsk_is) );
		TRYCXX( ISDestroy(&integralsk_is) );

	} /* endfor through clusters */


	this->timer_solve.stop();

	LOG_FUNC_END
}

/* ----------------------- external content */
void EntropySolverNewton<PetscVector>::ExternalContent::compute_gradient(Vec &g_Vec, Vec &integrals_Vec, Vec &moments_Vec) {
	LOG_FUNC_BEGIN

	int n;
	TRYCXX( VecGetSize(moments_Vec, &n) );

    double *g_arr;
    double *integrals_arr;
    double *moments_arr;

	TRYCXX( VecGetArray(g_Vec, &g_arr) );
	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );

    /* compute gradient */

	for (int j = 0; j < n; j++){
        g_arr[j] = moments_arr[j] - integrals_arr[1+j]/integrals_arr[0];
	}

	TRYCXX( VecRestoreArray(g_Vec, &g_arr) );
	TRYCXX( VecRestoreArray(integrals_Vec, &integrals_arr) );
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	LOG_FUNC_END
}

void EntropySolverNewton<PetscVector>::ExternalContent::compute_hessian(Vec &integrals_Vec) {
	LOG_FUNC_BEGIN

	int n;
	TRYCXX( MatGetSize(H_petsc, &n, NULL) );

    double *integrals_arr;

	TRYCXX( VecGetArray(integrals_Vec, &integrals_arr) );

    /* fill Hessian matrix */
	int counter = 1+n;
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			double temp = integrals_arr[counter]/integrals_arr[0] - integrals_arr[1+i]*integrals_arr[1+j]/(integrals_arr[0]*integrals_arr[0]);
			counter++;
			TRYCXX( MatSetValue(H_petsc, i, j, temp, INSERT_VALUES) );
			TRYCXX( MatSetValue(H_petsc, j, i, temp, INSERT_VALUES) );
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

//#endif
