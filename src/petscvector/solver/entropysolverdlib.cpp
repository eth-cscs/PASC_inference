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
	externalcontent = new ExternalContent(this->integration_eps, this->integration_type, this->debug_print_content);

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

	externalcontent->Fs = new double[entropydata->get_K()];

	LOG_FUNC_END
}


template<> 
EntropySolverDlib<PetscVector>::~EntropySolverDlib(){
	LOG_FUNC_BEGIN

	/* destroy auxiliary vectors of powers */
	free(externalcontent->x_powers_Vecs);
	//TODO: VecDestroy to x_powers_Vecs?

	free(externalcontent->Fs);

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
	int xdim = entropydata->get_xdim();

	/* Anna knows the purpose of this number */
	double eps = 0.0;

	/* prepare objects for Dlib */
	int nmb_of_moments = entropydata->get_number_of_moments();
    column_vector Mom(nmb_of_moments-1);
    column_vector starting_point(nmb_of_moments-1);

	/* stuff for PETSc to Dlib */
	Vec moments_Vec = moments->get_vector();
	Vec lambda_Vec = entropydata->get_lambda()->get_vector();
	double *moments_arr;
	double *lambda_arr;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
	TRYCXX( VecGetArray(lambda_Vec, &lambda_arr) );

	/* mom_powers: from D_matrix to dlib matrix, without first zero vector */
    dlib::matrix<double> mom_powers;
    mom_powers.set_size(nmb_of_moments-1, xdim);

	Vec matrix_D_Vec = entropydata->get_matrix_D()->get_vector();
	double *matrix_D_arr;
	TRYCXX( VecGetArray(matrix_D_Vec, &matrix_D_arr) );
	for(int D_row_idx=0; D_row_idx < nmb_of_moments-1; D_row_idx++){
		for(int D_col_idx=0; D_col_idx < xdim; D_col_idx++){
			mom_powers(D_row_idx,D_col_idx) = matrix_D_arr[(D_row_idx+1)*xdim + D_col_idx];
		}
	}
	TRYCXX( VecRestoreArray(matrix_D_Vec, &matrix_D_arr) );

	/* through all clusters */
	for(int k = 0; k < K; k++){
		/* Mom: from PETSc vector to Dlib column_vector, without first component */
		for(int idx=0;idx < nmb_of_moments-1;idx++){
			Mom(idx) = moments_arr[k*nmb_of_moments + idx+1];
		}

		/* prepare lambda-functions for Dlib */
		auto get_functions_obj_lambda = [&](const column_vector& x)->double { return externalcontent->get_functions_obj(x, Mom, eps, Km, mom_powers);};
		auto get_functions_grad_lambda = [&](const column_vector& x)->column_vector { return externalcontent->get_functions_grad(x, Mom, Km, mom_powers);};
		auto get_functions_hess_lambda = [&](const column_vector& x)->dlib::matrix<double> { return externalcontent->get_functions_hess(x, Mom, Km, mom_powers);};

		/* initial value forms starting_point, take solution from previous iteration */
//		for(int idx=0;idx<nmb_of_moments-1;idx++){
//			starting_point(idx) = lambda_arr[k*nmb_of_moments+idx+1];
//		}
		starting_point = 0.0;

		/* print cluster info */
		if(debug_print_it){
			coutMaster << "cluster = " << k << std::endl;
		}
		
		if(sum(Mom) != 0){
			/* solve using Dlib magic */
			if(debug_print_it){
				coutMaster << "- start to solve the problem using Dlib solver" << std::endl;
				/* give iteration info */
/*
				dlib::find_min_box_constrained(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
*/

/*
				dlib::find_min(dlib::bfgs_search_strategy(),
                             dlib::objective_delta_stop_strategy(this->eps).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point,-1e12);
*/
/*
				dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                             dlib::objective_delta_stop_strategy(this->eps).be_verbose(),
                             get_functions_obj_lambda, starting_point,-1e12);
*/

/*
				dlib::find_min(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point,-1e12);
*/

				dlib::find_min(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps,this->maxit).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point,-1e12);

				coutMaster << "- solver finished" << std::endl;
			} else {
				/* be quite */
				dlib::find_min_box_constrained(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(this->eps),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
			
			}
		} else {
			coutMaster << "ERROR: sum(moments) = 0" << std::endl;
			
			for(int idx=0;idx<nmb_of_moments-1;idx++){
				starting_point(idx) = 1000000;//std::numeric_limits<double>::max();//lambda_arr[k*nmb_of_moments+idx+1];
			}
		}
		
		externalcontent->Fs[k] = externalcontent->get_F();

		/* print cluster info */
/*		if(debug_print_it){
			coutMaster << "solution: " << std::endl;
			coutMaster << starting_point << std::endl;
			coutMaster << "----------------------" << std::endl;
		}
*/
		/* store lambda (solution): from Dlib to Petsc */
		lambda_arr[k*nmb_of_moments] = 1.0;
		for(int idx=0;idx<nmb_of_moments-1;idx++){
			lambda_arr[k*nmb_of_moments+idx+1] = starting_point(idx);
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
	int number_of_moments = entropydata->get_number_of_moments();

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
	int xdim = entropydata->get_xdim();
	int number_of_moments = entropydata->get_number_of_moments();

	/* update gamma_solver data - prepare new linear term */
	/* theta includes all moments */
	double *lambda_arr;
	TRYCXX( VecGetArray(entropydata->get_lambda()->get_vector(), &lambda_arr) );
	
	/* mom_powers - exponents */
	Vec matrix_D_Vec = entropydata->get_matrix_D()->get_vector();
	double *matrix_D_arr;
	TRYCXX( VecGetArray(matrix_D_Vec, &matrix_D_arr) );	

	/* index set of the data for given dimension component 1,...,xdim */
	IS xn_is;

	/* indexes of appropriate components in residuum 1,...,K */
	IS gammak_is;
	Vec gammak_Vec;

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&temp_Vec) );
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
		logF = log(externalcontent->Fs[k]);

		this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(residuum_Vec, gammak_is, &gammak_Vec) );

		/* multiply with correspoinding computed lagrange multiplier and add it to residuum */
		TRYCXX( VecSet(gammak_Vec, logF) );
	
		TRYCXX( VecRestoreSubVector(residuum_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	int D_value;
	for(int D_row_idx=0; D_row_idx < number_of_moments; D_row_idx++){ /* go through all rows of matrix D */

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
	TRYCXX( VecRestoreArray(matrix_D_Vec, &matrix_D_arr) );	
	TRYCXX( VecRestoreArray(entropydata->get_lambda()->get_vector(), &lambda_arr) );
		
	LOG_FUNC_END
}

template<>
double EntropySolverDlib<PetscVector>::get_integration_time() const {
	return externalcontent->integration_time;
}


/* -------------- External content -------------- */

EntropySolverDlib<PetscVector>::ExternalContent::ExternalContent(double new_integration_eps, int integration_type, bool debug_print_content) {
	this->integration_eps = new_integration_eps;
	this->integration_time = 0.0;
	this->integration_type = integration_type;
	this->debug_print_content = debug_print_content;
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

double EntropySolverDlib<PetscVector>::ExternalContent::get_functions_obj(const column_vector& _LM, const column_vector& _Mom, double eps, int k, const dlib::matrix<double>& mom_powers){
	column_vector LM = _LM;
	this->cLM = _LM;
	column_vector Mom = _Mom;
	dlib::matrix<double> D = mom_powers;
	
	/* number of variables */
	long n = D.nr();
    
//	column_vector I(n); /* theoretical moments */
	column_vector grad(n); /* gradient */
	dlib::matrix<double> hess; /* hessian */
	hess.set_size(n,n);
    
	int number_of_integrals = 1 + n + (int)(0.5*n*(n+1));

	Integrator integrator(this->integration_type, D.nc(), number_of_integrals);

	/* setting to compute normalization constant */
	ExtraParameters xp(k, Mom, LM, 0.0, 0.0, D, 0, 0);
	integrator.USERDATA = &xp;

	cubareal *computed_integrals;
	computed_integrals = integrator.compute();
    
	double F_ = computed_integrals[0];

	if(isinf(F_)){
		coutMaster << "ERROR: infinite objective function, F=" << F_ << std::endl;
	}
    
	/* modify settings to compute theoretical moments */
	for (int j = 0; j < n; j++){
		grad(j) = Mom(j) - computed_integrals[1+j]/F_;
	}

	/* modify settings to compute extra theoretical moments for hess */
	double temp = 0.0;
	int counter = 1+n;
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			temp = computed_integrals[counter]/F_ - computed_integrals[1+i]*computed_integrals[1+j]/(F_*F_);
			counter++;

			hess(i,j) = temp;
			hess(j,i) = temp;
		}
	}

	this->cgrad = grad;
	this->chess = hess;
	this->cF = F_;
	this->integration_time += integrator.get_time();

	if(debug_print_content){
		std::cout << "_F:" << std::endl;
		std::cout << F_ << std::endl; 
		std::cout << std::endl;

		std::cout << "Mom" << std::endl;
		std::cout << Mom << std::endl; 

		std::cout << "objective" << std::endl;
		std::cout << trans(Mom)*LM + log(F_) << std::endl; 

		std::cout << "LM" << std::endl;
		std::cout << _LM << std::endl;

		std::cout << "gradient" << std::endl;
		std::cout << this->cgrad << std::endl;

		std::cout << "Hessian matrix" << std::endl;
		std::cout << this->chess << std::endl;
	}

	return trans(Mom)*LM + log(F_);// + eps*sum(LM);
}

column_vector EntropySolverDlib<PetscVector>::ExternalContent::get_functions_grad(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers){
    if (_LM != cLM){
		coutMaster << "ERROR: gradient - lagrange multipliers changed" << std::endl;
	}
    return this->cgrad;
}

dlib::matrix<double> EntropySolverDlib<PetscVector>::ExternalContent::get_functions_hess(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers){
    if (_LM != cLM){
		coutMaster << "ERROR: hessian - lagrange multipliers changed" << std::endl;
	}
	
    return this->chess;
}


double EntropySolverDlib<PetscVector>::ExternalContent::get_F() const{
	return this->cF;
}

template<> EntropySolverDlib<PetscVector>::ExternalContent * EntropySolverDlib<PetscVector>::get_externalcontent() const {
	return externalcontent;
}


/* ------------ ExtraParameters ----------------- */
EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters(){
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

EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters(int _k, column_vector _Mom, column_vector _LM, double _L0, double _eps, dlib::matrix<double> _D, int _type, int _order){
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

void EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::Copy(ExtraParameters &_ExtraParameters){
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

EntropySolverDlib<PetscVector>::ExternalContent::ExtraParameters::~ExtraParameters(){
}


/* ------------ Integrator ----------------- */
int EntropySolverDlib<PetscVector>::ExternalContent::Integrator::Integrand(const int *ndim, const cubareal xx[],
                     const int *ncomp, cubareal ff2[], void *userdata) {
    

    //double p1 = xx[0];
    //double p2 = xx[1];
    //double p3 = xx[2];

    ExtraParameters* xp = (ExtraParameters*)userdata;
    dlib::matrix<double> D = xp->D;
    long d = D.nc();
    long n = D.nr();

//    int type = xp->type;
    column_vector LM = xp->LM;

    double V = 0.0;
    double p = 0.0;
    
    for (int i = 0; i < n; i++){
        p = 1.0;
        for (int j = 0; j < d; j++){
            p = p*pow(xx[j], D(i,j));
		}
        V = V - p*LM(i);
    }

	/* ff2[0] - type = 0 */
    ff2[0] = exp(V);
    
    /* ff2[1-n] - for gradient, type =2 */
	for(int order = 0; order < n; order++){
		p = 1.0;
		for (int j = 0; j < d; j++){
			p = p*pow(xx[j], D(order,j));
		}
        ff2[1+order] = p*exp(V);
    }

	/* ff2[n+1 - n+1+n*(n+1)/2] - for Hessian, type = 3 */
	int counter = 1+n;
	for(int order = 0; order < n; order++){
		for(int order2 = order; order2 < n; order2++){
			p = 1.0;
			row_vector t = rowm(D,order) + rowm(D,order2);
			for (int j = 0; j < d; j++){
				p = p*pow(xx[j], t(j));
			}
			ff2[counter] = p*exp(V);
			counter++;
		}
	}
	
    return 0;
}

EntropySolverDlib<PetscVector>::ExternalContent::Integrator::Integrator(int integration_type, int ndim, int ncomp)
{
    //all this paramterers are from example file demo-c.c
    NDIM = ndim; 
    NCOMP = ncomp;
    NVEC = 1;
    EPSREL = 1e-12;//1e-3;
    EPSABS = 1e-12;
    VERBOSE = 0; //log output
    LAST = 4;
    SEED = 0;
    MINEVAL = 0;
    MAXEVAL = 50000;
    NSTART = 1000;
    NINCREASE = 500;
    NBATCH = 1000;
    GRIDNO = 0;
    STATEFILE = NULL;
    SPIN = NULL;
    NNEW = 1000;
    NMIN = 2;
    FLATNESS = 25.0;
    USERDATA = NULL; //this is to pass extra parameters to integral
    
    KEY1 = 47;
    KEY2 = 1;
    KEY3 = 1;
    MAXPASS = 5;
    BORDER = 0.0;
    MAXCHISQ = 10.0;
    MINDEVIATION = 0.25;
    NGIVEN = 0;
    LDXGIVEN = NDIM;
    NEXTRA = 0;
    
    KEY = 0;
    
    this->integration_type = integration_type;
    this->integral = new cubareal[ncomp];
    this->error = new cubareal[ncomp];
    this->prob = new cubareal[ncomp];
    
    /* restart timer */
    timer.restart();
}

double EntropySolverDlib<PetscVector>::ExternalContent::Integrator::get_time() const {
	return timer.get_value_sum();
}

cubareal* EntropySolverDlib<PetscVector>::ExternalContent::Integrator::compute() {
	switch(this->integration_type){
		case 0: this->computeVegas(); break;
		case 1: this->computeSuave(); break;
		case 2: this->computeDivonne(); break;
		case 3: this->computeCuhre(); break;
	}
    return integral;
}

void EntropySolverDlib<PetscVector>::ExternalContent::Integrator::computeVegas(){
    //printf("-------------------- Vegas test --------------------\n");

	timer.start();
    
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
          EPSREL, EPSABS, VERBOSE, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, SPIN,
          &neval, &fail, integral, error, prob);

	timer.stop();
    
    //printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    //       neval, fail);
    //for( comp = 0; comp < NCOMP; ++comp )
    //    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //           (double)integral[comp], (double)error[comp], (double)prob[comp]);
}

void EntropySolverDlib<PetscVector>::ExternalContent::Integrator::computeSuave(){
//    printf("\n-------------------- Suave test --------------------\n");

	timer.start();
    
    Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST, SEED,
              MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);

	timer.stop();
    
//    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
//               nregions, neval, fail);
//    for( comp = 0; comp < NCOMP; ++comp )
//        printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
//               (double)integral[comp], (double)error[comp], (double)prob[comp]);
}

void EntropySolverDlib<PetscVector>::ExternalContent::Integrator::computeDivonne(){
//    printf("\n------------------- Divonne test -------------------\n");

	timer.start();
    
    Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                BORDER, MAXCHISQ, MINDEVIATION,
                NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

	timer.stop();
    
//    printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
//            nregions, neval, fail);
//    for( comp = 0; comp < NCOMP; ++comp )
//        printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
//              (double)integral[comp], (double)error[comp], (double)prob[comp]);
}

void EntropySolverDlib<PetscVector>::ExternalContent::Integrator::computeCuhre(){
//    printf("\n-------------------- Cuhre test --------------------\n");

	timer.start();
    
    Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST,
              MINEVAL, MAXEVAL, KEY,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);

	timer.stop();
    
//    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
//              nregions, neval, fail);
//    for( comp = 0; comp < NCOMP; ++comp )
//        printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
//               (double)integral[comp], (double)error[comp], (double)prob[comp]);
}

EntropySolverDlib<PetscVector>::ExternalContent::Integrator::~Integrator(){
	free(this->integral);
	free(this->error);
	free(this->prob);
}


}
} /* end namespace */

#endif /* Cuba */
#endif /* Dlib */
