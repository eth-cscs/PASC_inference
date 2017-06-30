#include "external/petscvector/solver/entropysolverdlib.h"

#ifdef USE_DLIB
#ifdef USE_CUBA
/* if DLib and Cuba are not used, then this class is quite useless */

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

	/* prepare and compute auxiliary array of powers */
	// TODO: aaaah! bottleneck, maybe can be performed in different way (I hope so)
	int Km = entropydata->get_Km();
	Vec x_Vec = entropydata->get_x()->get_vector();
	externalcontent->x_powers_Vecs = new Vec[Km+1]; /* 0 = x^0, 1 = x^1 ... Km = x^Km */
	TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[0])) );
	TRYCXX( VecSet(externalcontent->x_powers_Vecs[0], 1.0) );
	for(int km = 1; km <= Km;km++){
		TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[km])) );
		TRYCXX( VecPointwiseMult(externalcontent->x_powers_Vecs[km], externalcontent->x_powers_Vecs[km-1], x_Vec) );
	}

	LOG_FUNC_END
}


template<> 
EntropySolverDlib<PetscVector>::~EntropySolverDlib(){
	LOG_FUNC_BEGIN

	/* destroy auxiliary vectors of powers */
	free(externalcontent->x_powers_Vecs);
	//TODO: VecDestroy to x_powers_Vecs?

	/* destroy external content */
	free(externalcontent);

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
	int nmb_of_moments = entropydata->get_number_of_moments();
    column_vector Mom(nmb_of_moments);
    column_vector starting_point(nmb_of_moments);

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
		for(int idx=0;idx < nmb_of_moments;idx++){
			Mom(idx) = moments_arr[k*nmb_of_moments + idx];
		}

		/* prepare lambda-functions for Dlib */
		auto get_functions_obj_lambda = [&](const column_vector& x)->double { return externalcontent->get_functions_obj(x, Mom, eps);};
		auto get_functions_grad_lambda = [&](const column_vector& x)->column_vector { return externalcontent->get_functions_grad(x, Mom, Km);};
		auto get_functions_hess_lambda = [&](const column_vector& x)->dlib::matrix<double> { return externalcontent->get_functions_hess(x, Mom, Km);};

		/* initial value forms starting_point */
		starting_point = 0.0;

		/* print cluster info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}
		
		/* solve using Dlib magic */
		if(debug_print_it){
			/* give iteration info */
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
		for(int idx=0;idx<nmb_of_moments;idx++){
			lambda_arr[k*nmb_of_moments+idx] = starting_point(idx);
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

	/* I assume that externalcontent->x_powers_Vecs is computed and constant */
	/* I assume that D_matrix is computed and prepared */

	int xdim = entropydata->get_xdim();

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments->get_vector();

	Vec gammak_Vec;
	IS gammak_is;

	Vec matrix_D_Vec = entropydata->get_matrix_D()->get_vector();

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, this->entropydata->get_decomposition()->get_Tlocal(), this->entropydata->get_decomposition()->get_T()) );
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

	double *matrix_D_arr;
	TRYCXX( VecGetArray(matrix_D_Vec, &matrix_D_arr) );

	double mysum, gammaksum;

	int D_value;

	for(int D_row_idx=0; D_row_idx < entropydata->get_number_of_moments(); D_row_idx++){ /* go through all rows of matrix D */

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
				moments_arr[k*entropydata->get_number_of_moments() + D_row_idx] = mysum/gammaksum;
			} else {
				moments_arr[k*entropydata->get_number_of_moments() + D_row_idx] = 0.0;
			}
	
			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );
		}
		
	}
	
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );
	TRYCXX( VecRestoreArray(matrix_D_Vec, &matrix_D_arr) );

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


/* ------------ ExtraParameters ----------------- */
EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters()
{
    k = 0;
    D = 0.0;
    eps= 0.0;
    Mom = 0.0;
    type = 0;
    LM = 0.0;
    L0 = 0.0;
    order = 0;
    order2 = 0;
}

EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters(int _k, column_vector _Mom, column_vector _LM, double _L0, double _eps, dlib::matrix<double> _D, int _type, int _order)
{
    k = _k;
    Mom = _Mom;
    LM = _LM;
    L0 = _L0;
    eps = _eps;
    D = _D;
    type = _type;
    order = _order;
    order2 = _order;
}

void EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::Copy(ExtraParameters &_ExtraParameters)
{
    k = _ExtraParameters.k;
    Mom = _ExtraParameters.Mom;
    eps = _ExtraParameters.eps;
    D = _ExtraParameters.D;
    type = _ExtraParameters.type;
    LM = _ExtraParameters.LM;
    L0 = _ExtraParameters.L0;
    order = _ExtraParameters.order;
    order2 = _ExtraParameters.order2;
}

EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::~ExtraParameters()
{
}

}
} /* end namespace */

#endif /* Cuba */
#endif /* Dlib */
