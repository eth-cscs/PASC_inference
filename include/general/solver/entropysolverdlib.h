/** @file entropysolverdlib.h
 *  @brief Solver which solves problem with integrals using DLib library
 *
 *  @author Anna Marchenko & Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERDLIB_H
#define	PASC_ENTROPYSOLVERDLIB_H

#include "general/solver/generalsolver.h"
#include "general/data/entropydata.h"

/* include Dlib stuff */
#ifdef USE_DLIB
	#include "dlib/matrix.h"
	#include "dlib/numeric_constants.h"
	#include "dlib/numerical_integration.h"
	#include "dlib/optimization.h"

	/* Dlib column vector */
	typedef dlib::matrix<double,0,1> column_vector;
#endif

#define STOP_TOLERANCE 1e-06

namespace pascinference {
namespace solver {

/** \class EntropySolverDlib
 *  \brief Solver which solves problem with integrals using DLib library
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolverDlib: public GeneralSolver {
	protected:
		Timer timer_solve; /**< total solution time */
		Timer timer_compute_moments; /**< time of moment computation */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments */
		GeneralVector<VectorBase> *x_power; /**< temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< temp vector for storing power of x * gamma_k */

		#ifdef USE_DLIB
			/* functions for Dlib */
			static double gg(double y, int order, const column_vector& LM);
			double get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps);
			column_vector get_functions_grad(const column_vector& LM, const column_vector& Mom, int k);
			dlib::matrix<double> get_functions_hess(const column_vector& LM, const column_vector& Mom, int k);
		#endif

	public:

		EntropySolverDlib();
		EntropySolverDlib(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolverDlib();

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
namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
EntropySolverDlib<VectorBase>::EntropySolverDlib(){
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
EntropySolverDlib<VectorBase>::EntropySolverDlib(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	//TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolverDlib<VectorBase>::~EntropySolverDlib(){
	LOG_FUNC_BEGIN

	free(x_power);
	free(x_power_gammak);
	free(moments);

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolverDlib<VectorBase>::print(ConsoleOutput &output) const {
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
void EntropySolverDlib<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
void EntropySolverDlib<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	
	//TODO
	
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void EntropySolverDlib<VectorBase>::printcontent(ConsoleOutput &output) const {
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
void EntropySolverDlib<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve   = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_moments = " << this->timer_compute_moments.get_value_sum() << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolverDlib<VectorBase>::get_name() const {
	std::string return_value = "EntropySolverDlib<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
EntropyData<VectorBase> *EntropySolverDlib<VectorBase>::get_data() const {
	return entropydata;
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::compute_moments() {
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::compute_residuum(GeneralVector<VectorBase> *residuum) const {
	LOG_FUNC_BEGIN

	//TODO
		
	LOG_FUNC_END
}


#ifdef USE_DLIB

template<class VectorBase>
double EntropySolverDlib<VectorBase>::gg(double y, int order, const column_vector& LM){
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    
    return pow(y,order)*(exp(-trans(LM)*z));
}

template<class VectorBase>
double EntropySolverDlib<VectorBase>::get_functions_obj(const column_vector& LM, const column_vector& Mom, double eps){
    /* compute normalization */
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    return dlib::trans(Mom)*LM + log(F_);// + eps*sum(LM);	
}

template<class VectorBase>
column_vector EntropySolverDlib<VectorBase>::get_functions_grad(const column_vector& LM, const column_vector& Mom, int k){
    column_vector grad(k);
    column_vector I(k);
    
    /* compute normalization */
    column_vector LMVec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, LMVec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    /* theoretical moments */
    int i = 0;
    while (i < k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, LMVec);};
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
dlib::matrix<double> EntropySolverDlib<VectorBase>::get_functions_hess(const column_vector& LM, const column_vector& Mom, int k){
    dlib::matrix<double> hess(k, k);
    
    column_vector I(2*k);
    
    //compute normalization
    column_vector LMVec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, LMVec);};//std::bind(gg, _1,  1, 2);
    double F_ = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    //theoretical moments
    int i = 0;
    while (i < 2*k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, LMVec);};
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

#endif /* USE_DLIB */

}
} /* end namespace */

#endif
