/** @file entropysolver.h
 *  @brief Solver which solves problem with integrals from Anna
 *
 *  @author Anna Marchenko & Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVER_H
#define	PASC_ENTROPYSOLVER_H

#include "pascinference.h"
#include "data/entropydata.h"

/* this code is for Dlib */
#ifndef USE_DLIB
 #error 'ENTROPYSOLVER is for DLIB'
#endif

/* include Dlib stuff */
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

#define STOP_TOLERANCE 1e-06

/* Dlib column vector */
typedef dlib::matrix<double,0,1> column_vector;

namespace pascinference {
namespace solver {

/** \class EntropySolver
 *  \brief Solver which solves problem with integrals from Anna
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolver: public GeneralSolver {
	protected:
		Timer timer_solve; /**< total solution time */
		Timer timer_compute_moments; /**< time of moment computation */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments */
		GeneralVector<VectorBase> *x_power; /**< temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< temp vector for storing power of x * gamma_k */

		/* functions for Dlib */
		double gg(double y, int order, column_vector& LM);
		double get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps);
		column_vector get_functions_grad(const column_vector& LM, const column_vector& Mom, int k);
		dlib::matrix<double> get_functions_hess(const column_vector& LM, const column_vector& Mom, int k);

	public:

		EntropySolver();
		EntropySolver(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printcontent(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		EntropyData<VectorBase> *get_data() const;

		void compute_moments();
		
		void compute_residuum(GeneralVector<VectorBase> *residuum) const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
EntropySolver<VectorBase>::EntropySolver(){
	LOG_FUNC_BEGIN
	
	entropydata = NULL;
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	LOG_FUNC_END
}

template<class VectorBase>
EntropySolver<VectorBase>::EntropySolver(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	

	/* prepare auxiliary vectors */
	x_power = new GeneralVector<PetscVector>(*entropydata->get_x());
	x_power_gammak = new GeneralVector<PetscVector>(*entropydata->get_x());
	
	/* create aux vector for the computation of moments */
	Vec moments_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&moments_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(moments_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(moments_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*entropydata->get_Km(),PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	this->moments = new GeneralVector<PetscVector>(moments_Vec);

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolver<VectorBase>::~EntropySolver(){
	LOG_FUNC_BEGIN

	free(x_power);
	free(x_power_gammak);
	free(moments);

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		coutMaster.push();
		entropydata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void EntropySolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << std::endl;
	output_global <<  " - eps:        " << this->eps << std::endl;
	output_global <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		output_global << "- data:" << std::endl;
		coutMaster.push();
		entropydata->print(output_global,output_local);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	
	//TODO
	
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void EntropySolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(entropydata){
		output << "- data:" << std::endl;
		coutMaster.push();
		entropydata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve   = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_moments = " << this->timer_compute_moments.get_value_sum() << std::endl;


	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolver<VectorBase>::get_name() const {
	std::string return_value = "EntropySolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
EntropyData<VectorBase> *EntropySolver<VectorBase>::get_data() const {
	return entropydata;
}

template<class VectorBase>
void EntropySolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->compute_moments();
	
//	coutMaster << "Moments: " << *moments << std::endl;

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

	/* prepare lambda-functions for Dlib */
	auto get_functions_obj_lambda = [&](const column_vector& x)->double { return get_functions_obj(x, Mom, eps);};
	auto get_functions_grad_lambda = [&](const column_vector& x)->column_vector { return get_functions_grad(x, Mom, Km);};
	auto get_functions_hess_lambda = [&](const column_vector& x)->dlib::matrix<double> { return get_functions_hess(x, Mom, Km);};

	/* through all clusters */
	for(int k = 0; k < K; k++){
		/* Mom: from PETSc vector to Dlib column_vector */
		for(int km=0;km<Km;km++){
			Mom(km) = moments_arr[k*Km+km];
		}

		/* initial value form starting_point */
		starting_point = 0.0;

		coutMaster << "k=" << k << ": running Dlib miracle" << std::endl;

		/* solve using Dlib magic */
		dlib::find_min_box_constrained(dlib::newton_search_strategy(get_functions_hess_lambda),
                             dlib::objective_delta_stop_strategy(STOP_TOLERANCE).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );

//		coutMaster << "something computed" << std::endl;

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
void EntropySolver<PetscVector>::compute_moments() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start(); 

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec x_power_Vec = x_power->get_vector();
	Vec x_power_gammak_Vec = x_power_gammak->get_vector();
	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments->get_vector();
	
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

	this->timer_compute_moments.stop(); 
	
	LOG_FUNC_END
}


template<class VectorBase>
double EntropySolver<VectorBase>::gg(double y, int order, column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    
    return pow(y,order)*(exp(-trans(LM)*z));
}

template<class VectorBase>
double EntropySolver<VectorBase>::get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps){
    /* compute normalization */
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    return dlib::trans(Mom)*LM + log(F_);// + eps*sum(LM);	
}

template<class VectorBase>
column_vector EntropySolver<VectorBase>::get_functions_grad(const column_vector& LM, const column_vector& Mom, int k){
    column_vector grad(k);
    column_vector I(k);
    
    /* compute normalization */
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    /* theoretical moments */
    int i = 0;
    while (i < k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, Vec);};
        I(i) = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
        i++;
    }
    
    for (int i = 0; i < k; ++i)
        grad(i) = Mom(i) - I(i)/F_;
    
//    double L1 = grad(0);
//    double L2 = grad(1);
    return grad;
}

template<class VectorBase>
dlib::matrix<double> EntropySolver<VectorBase>::get_functions_hess(const column_vector& LM, const column_vector& Mom, int k){
    dlib::matrix<double> hess(k, k);
    
    column_vector I(2*k);
    
    //compute normalization
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    //theoretical moments
    int i = 0;
    while (i < 2*k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, Vec);};
        I(i) = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
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

template<>
void EntropySolver<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum) const {
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

}
} /* end namespace */

#endif
