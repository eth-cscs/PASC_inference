#include "external/petscvector/solver/entropysolverdlib.h"

#ifdef USE_DLIB
/* if DLib is not used, then this class is quite useless */

namespace pascinference {
namespace solver {

template<>
EntropySolverDlib<PetscVector>::EntropySolverDlib(EntropyData<PetscVector> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->integration_eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	/* prepare auxiliary vectors */
//	x_power = new GeneralVector<PetscVector>(*entropydata->get_x());
//	x_power_gammak = new GeneralVector<PetscVector>(*entropydata->get_x());
	
	/* create aux vector for the computation of moments */
	Vec moments_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&moments_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(moments_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(moments_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(moments_Vec,get_K()*get_number_of_moments(),PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	this->moments = new GeneralVector<PetscVector>(moments_Vec);

	/* prepare external content with PETSc-DLIB stuff */
	externalcontent = new ExternalContent(this->integration_eps);

	LOG_FUNC_END
}

template<>
void EntropySolverDlib<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	this->compute_moments();

	if(debug_print_moments){
		coutMaster << "Moments: " << *moments << std::endl;
	}

	this->timer_solve.start(); 

	/* get dimensions */
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* Anna knows the purpose of this number */
	double eps = 0.0;

	/* prepare objects for Dlib */
    column_vector Mom(Km);
    column_vector starting_point(Km);

	/* stuff for PETSc to Dlib */
	Vec moments_Vec = moments->get_vector();
	Vec lambda_Vec = entropydata->get_lambda()->get_vector();
	double *moments_arr;
	double *lambda_arr;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
	TRYCXX( VecGetArray(lambda_Vec, &lambda_arr) );

	/* through all clusters */
	for(int k = 0; k < K; k++){
		/* Mom: from PETSc vector to Dlib column_vector */
		for(int km=0;km<Km;km++){
			Mom(km) = moments_arr[k*Km+km];
		}

		/* prepare lambda-functions for Dlib */
		auto get_functions_obj_lambda = [&](const column_vector& x)->double { return externalcontent->get_functions_obj(x, Mom, eps);};
		auto get_functions_grad_lambda = [&](const column_vector& x)->column_vector { return externalcontent->get_functions_grad(x, Mom, Km);};
		auto get_functions_hess_lambda = [&](const column_vector& x)->dlib::matrix<double> { return externalcontent->get_functions_hess(x, Mom, Km);};

		/* initial value form starting_point */
		starting_point = 0.0;

		/* print cluster info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}
		
		/* solve using Dlib magic */
		if(debug_print_it){
			/* give it info */
			dlib::find_min_box_constrained(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
		} else {
			/* be quite */
			dlib::find_min_box_constrained(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
			
		}

		/* store lambda (solution): from Dlib to Petsc */
		for(int km=0;km<Km;km++){
			lambda_arr[k*Km+km] = starting_point(km);
		}
		
	} /* endfor through clusters */

	TRYCXX( VecRestoreArray(lambda_Vec, &lambda_arr) );
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );
	
	this->timer_solve.stop(); 

	LOG_FUNC_END
}

template<>
void EntropySolverDlib<PetscVector>::compute_moments() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start(); 

	Vec x_Vec = entropydata->get_x()->get_vector();

	/* prepare and copute auxiliary array of powers */
	// TODO: aaaah!
//	Vec *x_powers_Vecs = new Vec(get_Km()+1); /* 0 = x^0, 1 = x^1 ... Km = x^Km */
//	TRYCXX( VecDuplicate(x_Vec, &(x_powers_Vecs[0])) );
//	TRYCXX( VecSet(x_powers_Vecs[0], 1.0) );
//	for(km = 1; km <= get_Km();km++){
//		TRYCXX( VecPointwiseMult(x_powers_Vecs[km], x_powers_Vecs[km-1], x_Vec) );
//	}

	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments->get_vector();

	Vec gammak_Vec;
	IS gammak_is;

	/* temp = x_n^D */
	/* temp2 = (temp*temp*...) */
	/* temp = gamma_k*temp2 */
	/* mom = sum(temp) */
	Vec temp_Vec;
	Vec temp2_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, this->entropydata->get_decomposition()->get_Tlocal(), this->entropydata->get_decomposition()->get_T()) );
	TRYCXX( VecSetFromOptions(temp_Vec) );	
	TRYCXX( VecDuplicate(temp_Vec, &temp2_Vec));
	
	double *moments_arr, mysum, gammaksum;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );

//	for(int idx=0; idx < entropydata->get_number_of_moments(); idx++){
		
//		for(int k=0;k<entropydata->get_K();k++){
			/* get gammak */
//			this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
//			TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

			/* compute x_power_gammak */
//			TRYCXX( VecPointwiseMult(x_power_gammak_Vec, gammak_Vec, x_power_Vec) ); /* x_power_gammak = x_power.*gammak */

			/* compute gammaksum */
//			TRYCXX( VecSum(gammak_Vec, &gammaksum) );
//			TRYCXX( VecSum(x_power_gammak_Vec, &mysum) );

			/* store computed moment */
//			if(gammaksum != 0){
//				moments_arr[k*this->entropydata->get_Km() + km] = mysum/gammaksum;
//			} else {
//				moments_arr[k*this->entropydata->get_Km() + km] = 0.0;
//			}
	
//			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
//			TRYCXX( ISDestroy(&gammak_is) );	
//		}
		
//		TRYCXX( VecPointwiseMult(x_power_Vec, x_Vec, x_power_Vec) ); /* x_power = x_power.*x */
//	}
	
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	this->timer_compute_moments.stop(); 
	
	LOG_FUNC_END
}

template<>
void EntropySolverDlib<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
	LOG_FUNC_BEGIN

	int T = entropydata->get_T();
	int Tlocal = entropydata->get_decomposition()->get_Tlocal();
	int K = entropydata->get_K();
	int Km = entropydata->get_Km();

	/* lambda vector for Dlib integration */
	column_vector lambda_Dlib(Km);
    auto mom_function = [&](double x)->double { return externalcontent->gg(x, 0, lambda_Dlib);};
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
			lambda_Dlib(km) = lambda_arr[k*Km+km];
		}
		F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
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


/* External content */
EntropySolverDlib<PetscVector>::ExternalContent::ExternalContent(double new_integration_eps) {
	this->integration_eps = new_integration_eps;
}

double EntropySolverDlib<PetscVector>::ExternalContent::gg(double y, int order, const column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    
    return pow(y,order)*(exp(-trans(LM)*z));
}

double EntropySolverDlib<PetscVector>::ExternalContent::get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps){
    /* compute normalization */
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
    
    return dlib::trans(Mom)*LM + log(F_);// + eps*sum(LM);	
}

column_vector EntropySolverDlib<PetscVector>::ExternalContent::get_functions_grad(const column_vector& LM, const column_vector& Mom, int k){
    column_vector grad(k);
    column_vector I(k);
    
    /* compute normalization */
    column_vector LMVec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, LMVec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
    
    /* theoretical moments */
    int i = 0;
    while (i < k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, LMVec);};
        I(i) = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
        i++;
    }
    
    for (int i = 0; i < k; ++i)
        grad(i) = Mom(i) - I(i)/F_;
    
//    double L1 = grad(0);
//    double L2 = grad(1);
    return grad;
}

dlib::matrix<double> EntropySolverDlib<PetscVector>::ExternalContent::get_functions_hess(const column_vector& LM, const column_vector& Mom, int k){
    dlib::matrix<double> hess(k, k);
    
    column_vector I(2*k);
    
    //compute normalization
    column_vector LMVec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, LMVec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
    
    //theoretical moments
    int i = 0;
    while (i < 2*k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, LMVec);};
        I(i) = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->integration_eps);
        i++;
    }
    
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            hess(i,j) = I(i+j+1)/F_ - I(i)*I(j)/(F_*F_);
    
//    double L1 = hess(0,0);
//    double L2 = hess(0,1);
//    double L3 = hess(1,0);
//    double L4 = hess(1,1);
    return hess;
}

template<> EntropySolverDlib<PetscVector>::ExternalContent * EntropySolverDlib<PetscVector>::get_externalcontent() const {
	return externalcontent;
}


}
} /* end namespace */

#endif
