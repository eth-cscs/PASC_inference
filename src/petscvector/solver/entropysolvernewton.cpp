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

	/* prepare and compute auxiliary array of powers */
	// TODO: aaaah! bottleneck, maybe can be performed in different way (I hope so)
	int Km = entropydata->get_Km();
	Vec x_Vec = entropydata->get_x()->get_vector();
	externalcontent->x_powers_Vecs = new Vec[Km+1]; /* 0 = x^0, 1 = x^1 ... Km = x^Km */
	TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[0])) );
	TRYCXX( VecSet(externalcontent->x_powers_Vecs[0], 1.0) );
	TRYCXX( VecAssemblyBegin(externalcontent->x_powers_Vecs[0]));
	TRYCXX( VecAssemblyEnd(externalcontent->x_powers_Vecs[0]));

	for(int km = 1; km <= Km;km++){
		TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[km])) );
		TRYCXX( VecPointwiseMult(externalcontent->x_powers_Vecs[km], externalcontent->x_powers_Vecs[km-1], x_Vec) );
		TRYCXX( VecAssemblyBegin(externalcontent->x_powers_Vecs[km]));
		TRYCXX( VecAssemblyEnd(externalcontent->x_powers_Vecs[km]));
	}
	
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
	
	/* destroy auxiliary vectors of powers */
	free(externalcontent->x_powers_Vecs);
	//TODO: VecDestroy to x_powers_Vecs?

	/* destroy external content */
	free(externalcontent);

	LOG_FUNC_END
}

template<> 
void EntropySolverNewton<PetscVector>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN


	int n = this->number_of_moments-1;
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
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*n,PETSC_DECIDE) );
	TRYCXX( VecSetSizes(integrals_Vec,entropydata->get_K()*number_of_integrals,PETSC_DECIDE) );
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
	int number_of_moments = entropydata->get_number_of_moments();
	int n = number_of_moments-1;
	int number_of_integrals = 1 + n + (int)(0.5*n*(n+1));
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* get PETSc vecs - LOCAL */
	Vec moments_Vec = moments_data->get_vector();
	Vec x_Vec = entropydata->get_lambda()->get_vector(); /* x:=lambda is unknowns */
	Vec g_Vec = g->get_vector();
	Vec delta_Vec = delta->get_vector();
	Vec integrals_Vec = integrals->get_vector();

	/* stuff for clusters (subvectors) - LOCAL */
	IS x_is;
	IS moments_is;
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


	//TODO:temp
	int lambda_size;
	TRYCXX( VecGetSize(x_Vec, &lambda_size) );
	int moments_size;
	TRYCXX( VecGetSize(moments_Vec, &moments_size) );
	int integrals_size;
	TRYCXX( VecGetSize(integrals_Vec, &integrals_size) );
	std::cout << "x_size: " << lambda_size << std::endl;
	std::cout << "moments_size: " << moments_size << std::endl;
	std::cout << "integrals_size: " << integrals_size << std::endl;

	/* through all clusters */
	for(int k = 0; k < K; k++){
		
		/* print iteration info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}
		
		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, n, k*number_of_moments+1, 1, &x_is) );
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, n, k*n, 1, &moments_is) );
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, number_of_integrals, k*number_of_integrals, 1, &integralsk_is) ); 
	
		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(x_Vec, x_is, &xk_Vec) );
		TRYCXX( VecGetSubVector(moments_Vec, moments_is, &momentsk_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* -------------- Newton algorithm (for k-th problem) ------------ */
		it = 0;
		itAxb = 0;
		
		/* compute integrals and gradient */
		this->timer_integrate.start();
		 TRYCXX( VecGetArray(xk_Vec,&xk_arr));
		 TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
 		  entropyintegration->compute(integralsk_arr, xk_arr, number_of_integrals);
		 TRYCXX( VecRestoreArray(xk_Vec,&xk_arr));
		 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
		this->timer_integrate.stop();
		this->timer_g.start();
//		 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
		this->timer_g.stop();
		this->timer_fs.start();
//		 fx = externalcontent->compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
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
			  entropyintegration->compute(integralsk_arr, xk_arr, number_of_integrals);
			 TRYCXX( VecRestoreArray(xk_Vec,&xk_arr));
			 TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));
			this->timer_integrate.stop();
			this->timer_g.start();
//			 externalcontent->compute_gradient(g_Vec, integralsk_Vec, momentsk_Vec);
			this->timer_g.stop();
			this->timer_fs.start();
//			 fx = externalcontent->compute_function_value(xk_Vec, integralsk_Vec, momentsk_Vec);
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
		TRYCXX( VecRestoreSubVector(x_Vec, x_is, &xk_Vec) );
		TRYCXX( VecRestoreSubVector(moments_Vec, moments_is, &momentsk_Vec) );
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* destroy index set */
		TRYCXX( ISDestroy(&x_is) );
		TRYCXX( ISDestroy(&moments_is) );
		TRYCXX( ISDestroy(&integralsk_is) );	

	} /* endfor through clusters */


	this->timer_solve.stop(); 

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::compute_moments_data() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start();

	int Tlocal = this->entropydata->get_decomposition()->get_Tlocal();
	int Rlocal = this->entropydata->get_decomposition()->get_Rlocal();

	int T = this->entropydata->get_decomposition()->get_T();
	int R = this->entropydata->get_decomposition()->get_R();

	/* I assume that externalcontent->x_powers_Vecs is computed and constant */
	/* I assume that D_matrix is computed and prepared */

	int xdim = entropydata->get_xdim();
	int number_of_moments = entropydata->get_number_of_moments();

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments_data->get_vector();

	Vec gammak_Vec;
	IS gammak_is;

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, Tlocal*Rlocal, T*R) );
	TRYCXX( VecSetFromOptions(temp_Vec) );

	/* temp2 = x_n^D */
	Vec temp2_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp2_Vec));

	/* temp3 = gammak*temp; */
	Vec temp3_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp3_Vec)); //TODO: I cannot reuse temp2 ?

	/* mom = sum(temp3) */
	IS xn_is;

	double *moments_arr;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );

	int *matrix_D_arr = entropydata->get_matrix_D();

	double mysum, gammaksum;

	int D_value;

	for(int D_row_idx=0; D_row_idx < number_of_moments; D_row_idx++){ /* go through all rows of matrix D */

		/* temp = 1 */
		TRYCXX( VecSet(temp_Vec, 1.0) );
		TRYCXX( VecAssemblyBegin(temp_Vec));
		TRYCXX( VecAssemblyEnd(temp_Vec));

		/* throught columns of D */
		for(int D_col_idx=0; D_col_idx < xdim; D_col_idx++){
			D_value = (int)matrix_D_arr[D_row_idx*xdim + D_col_idx];

			/* get x_n^D */
			this->entropydata->get_decomposition()->createIS_datan(&xn_is, D_col_idx);

			TRYCXX( VecGetSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );

			/* compute temp *= x_n^D */
			TRYCXX( VecPointwiseMult(temp_Vec, temp_Vec, temp2_Vec) );
			TRYCXX( VecAssemblyBegin(temp_Vec));
			TRYCXX( VecAssemblyEnd(temp_Vec));

			TRYCXX( VecRestoreSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );
			TRYCXX( ISDestroy(&xn_is) );
		}

		/* go throught clusters and multiply with coefficients */
		for(int k=0;k<entropydata->get_K();k++){
			/* get gammak */
			this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

			/* compute temp_Vec*gammak_Vec */
			TRYCXX( VecPointwiseMult(temp3_Vec, gammak_Vec, temp_Vec) ); /* x_power_gammak = x_power.*gammak */

			/* compute gammaksum */
			TRYCXX( VecSum(gammak_Vec, &gammaksum) );
			TRYCXX( VecSum(temp3_Vec, &mysum) );

			/* store computed moment */
			if(gammaksum != 0){
				moments_arr[k*number_of_moments + D_row_idx] = mysum/gammaksum;
			} else {
				coutMaster << "ERROR: norm(gammak) = 0" << std::endl;

				moments_arr[k*entropydata->get_number_of_moments() + D_row_idx] = 0.0;
			}

			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );
		}

	}

	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	this->timer_compute_moments.stop();

	LOG_FUNC_END
}

template<>
void EntropySolverNewton<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
	LOG_FUNC_BEGIN

	int T = entropydata->get_T();
	int Tlocal = entropydata->get_decomposition()->get_Tlocal();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();
	int xdim = entropydata->get_xdim();
	int number_of_moments = entropydata->get_number_of_moments();

	int n = number_of_moments-1;
	int number_of_integrals = 1 + n + (int)(0.5*n*(n+1));
	
	/* update gamma_solver data - prepare new linear term */
	/* theta includes all moments */
	double *lambda_arr;
	TRYCXX( VecGetArray(entropydata->get_lambda()->get_vector(), &lambda_arr) );

	double *integrals_arr;
	TRYCXX( VecGetArray(integrals->get_vector(), &integrals_arr) );

	/* mom_powers - exponents */
	int *matrix_D_arr = entropydata->get_matrix_D();

	/* index set of the data for given dimension component 1,...,xdim */
	IS xn_is;

	/* indexes of appropriate components in residuum 1,...,K */
	IS gammak_is;
	Vec gammak_Vec;

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, Tlocal, T) );
	TRYCXX( VecSetFromOptions(temp_Vec) );

	/* temp2 = x_n^D */
	Vec temp2_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp2_Vec));

	/* vector of residuum */
	Vec residuum_Vec = residuum->get_vector();
	TRYCXX( VecSet(residuum_Vec, 0.0) ); /* residuum = 0 */
	/* add log part to residuum */
	double logF;
	for(int k=0; k<K;k++){
		logF = log(integrals_arr[k*number_of_integrals]);

		this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(residuum_Vec, gammak_is, &gammak_Vec) );

		/* multiply with correspoinding computed lagrange multiplier and add it to residuum */
		TRYCXX( VecSet(gammak_Vec, logF) );

		TRYCXX( VecRestoreSubVector(residuum_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	int D_value;
	for(int D_row_idx=1; D_row_idx < number_of_moments; D_row_idx++){ /* go through all rows of matrix D */

		/* temp = 1 */
		TRYCXX( VecSet(temp_Vec, 1.0) );

		/* throught columns of D */
		for(int D_col_idx=0; D_col_idx < xdim; D_col_idx++){
			D_value = (int)matrix_D_arr[D_row_idx*xdim + D_col_idx];

			/* get x_n^D */
			this->entropydata->get_decomposition()->createIS_datan(&xn_is, D_col_idx);
			TRYCXX( VecGetSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );

			/* compute temp *= x_n^D */
			TRYCXX( VecPointwiseMult(temp_Vec, temp_Vec, temp2_Vec) );

			TRYCXX( VecRestoreSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );
			TRYCXX( ISDestroy(&xn_is) );
		}

		/* go throught clusters and multiply with lambda */
		for(int k=0;k<entropydata->get_K();k++){
			/* get gammak */
			this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(residuum_Vec, gammak_is, &gammak_Vec) );

			/* multiply with correspoinding computed lagrange multiplier and add it to residuum */
			TRYCXX( VecAXPY(gammak_Vec, lambda_arr[k*number_of_moments+D_row_idx], temp_Vec) );

			TRYCXX( VecRestoreSubVector(residuum_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );
		}

	}

	/* restore arrays */
	TRYCXX( VecRestoreArray(entropydata->get_lambda()->get_vector(), &lambda_arr) );
	TRYCXX( VecRestoreArray(integrals->get_vector(), &integrals_arr) );

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

	//TODO:temp
	int lambda_size;
	TRYCXX( VecGetSize(lambda_Vec, &lambda_size) );
	int moments_size;
	TRYCXX( VecGetSize(moments_Vec, &moments_size) );
	std::cout << "lambda_size: " << lambda_size << std::endl;
	std::cout << "moments_size: " << moments_size << std::endl;

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
