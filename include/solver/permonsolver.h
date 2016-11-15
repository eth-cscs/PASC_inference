/** @file permonsolver.h
 *  @brief solve QP with solvers implemented in Permon library
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_PERMONSOLVER_H
#define	PASC_PERMONSOLVER_H

#include <iostream>

#include "pascinference.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#ifndef USE_PETSCVECTOR
	#error 'BLOCKGRAPHSPARSEMATRIX is for PETSCVECTOR only, sorry'
#else
	typedef petscvector::PetscVector PetscVector;
#endif

#ifndef USE_PERMON
	#error 'PERMONSOLVER cannot be used without -DUSE_PERMON=ON'
#else
	#include "fllopqp.h" /* manipulation with quadratic programming problems (QP) */
	#include "fllopqps.h" /* manipulation with solvers (QPS) */
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif


namespace pascinference {
namespace solver {

/** \class PermonSolver
 *  \brief Interface with QP solvers implemented in Permon library
 *
 *  For solving QP on closed convex set described by separable simplexes.
*/
template<class VectorBase>
class PermonSolver: public QPSolver<VectorBase> {
	private:
		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();


		/** @brief allocate storage for auxiliary vectors used in computation
		* 
		*/
		void allocate_temp_vectors();

		/** @brief deallocate storage for auxiliary vectors used in computation
		* 
		*/
		void free_temp_vectors();

		/** @brief prepare permon QP and QPS from our qpdata
		* 
		*/
		void prepare_permon_objects();

		/* temporary vectors used during the solution process */
		GeneralVector<VectorBase> *Ad; 		/**< A*d */


		Timer timer_solve; 			/**< total solution time of used algorithm */

		QPData<VectorBase> *qpdata; /**< data on which the solver operates */
		double gP; 					/**< norm of projected gradient */

		QP qp;						/**< Quadratic Programming problem */
		QPS qps;					/**< Quadratic Programming solver */
		Mat BE;						/**< matrix of equality constraints */
		Vec cE;						/**< vector of equality constraints */
		
	public:
		/** @brief general constructor
		* 
		*/
		PermonSolver();

		/** @brief constructor based on provided data of problem
		* 
		* @param new_qpdata data of quadratic program
		*/
		PermonSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		* 
		*/
		~PermonSolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printtimer(ConsoleOutput &output) const;
		void printshort(std::ostringstream &header, std::ostringstream &values) const;
		void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		std::string get_name() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

template<class VectorBase>
void PermonSolver<VectorBase>::set_settings_from_console() {
/*	consoleArg.set_option_value("spgqpsolver_maxit", &this->maxit, SPGQPSOLVER_COEFF_DEFAULT_MAXIT);
	consoleArg.set_option_value("spgqpsolver_eps", &this->eps, SPGQPSOLVER_COEFF_DEFAULT_EPS);
	
	consoleArg.set_option_value("spgqpsolver_m", &this->m, SPGQPSOLVER_COEFF_DEFAULT_M);	
	consoleArg.set_option_value("spgqpsolver_gamma", &this->gamma, SPGQPSOLVER_COEFF_DEFAULT_GAMMA);	
	consoleArg.set_option_value("spgqpsolver_sigma1", &this->sigma1, SPGQPSOLVER_COEFF_DEFAULT_SIGMA1);	
	consoleArg.set_option_value("spgqpsolver_sigma2", &this->sigma2, SPGQPSOLVER_COEFF_DEFAULT_SIGMA2);	
	consoleArg.set_option_value("spgqpsolver_alphainit", &this->alphainit, SPGQPSOLVER_COEFF_DEFAULT_ALPHAINIT);	

	consoleArg.set_option_value("spgqpsolver_stop_normgp", &this->stop_normgp, SPGQPSOLVER_COEFF_STOP_NORMGP);
	consoleArg.set_option_value("spgqpsolver_stop_Anormgp", &this->stop_Anormgp, SPGQPSOLVER_COEFF_STOP_ANORMGP);
	consoleArg.set_option_value("spgqpsolver_stop_normgp_normb", &this->stop_normgp_normb, SPGQPSOLVER_COEFF_STOP_NORMGP_NORMB);
	consoleArg.set_option_value("spgqpsolver_stop_Anormgp_normb", &this->stop_Anormgp_normb, SPGQPSOLVER_COEFF_STOP_ANORMGP_NORMB);
	consoleArg.set_option_value("spgqpsolver_stop_difff", &this->stop_difff, SPGQPSOLVER_COEFF_STOP_DIFFF);	
	consoleArg.set_option_value("spgqpsolver_debugmode", &this->debugmode, SPGQPSOLVER_COEFF_DEFAULT_DEBUGMODE);
	
	debug_print_vectors = false;
	debug_print_scalars = false; 
	debug_print_it = false; 

	if(debugmode == 1){
		debug_print_it = true;
	}

	if(debugmode == 2){
		debug_print_it = true;
		debug_print_scalars = true;
	}

	consoleArg.set_option_value("tssolver_debug_print_it",		&debug_print_it, 		debug_print_it);
	consoleArg.set_option_value("tssolver_debug_print_vectors", &debug_print_vectors,	false);
	consoleArg.set_option_value("tssolver_debug_print_scalars", &debug_print_scalars, 	debug_print_scalars);
*/

}


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
PermonSolver<VectorBase>::PermonSolver(){
	LOG_FUNC_BEGIN

	qpdata = NULL;
	qp = NULL;
	qps = NULL;
	
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;
	
	this->fx = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	

	LOG_FUNC_END
}

template<class VectorBase>
PermonSolver<VectorBase>::PermonSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;
	
	/* allocate temp vectors */
	allocate_temp_vectors();
	
	this->it_sum = 0;
	this->hessmult_sum = 0;
	this->it_last = 0;
	this->hessmult_last = 0;

	this->fx = std::numeric_limits<double>::max();
	this->gP = std::numeric_limits<double>::max();

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	

	prepare_permon_objects();

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::prepare_permon_objects(){
	LOG_FUNC_BEGIN

	/* get local dimension of vector (all should have the same layout, so I take for instance b) */
//	GeneralVector<VectorBase> pattern = qpdata->get_b();
//	int global_size = pattern.size();
//	int local_size = pattern.local_size();

	/* assemble linear equality conditions */
/*	TRYCXX( MatCreate(PETSC_COMM_WORLD, &BE) );
	TRYCXX( MatSetSizes(BE,PETSC_DECIDE,,2,n);CHKERRV(ierr);
	ierr = MatSetFromOptions(B);CHKERRV(ierr);
	ierr = MatMPIAIJSetPreallocation(B,2,NULL,2,NULL);CHKERRV(ierr);
	ierr = MatSeqAIJSetPreallocation(B,2,NULL);CHKERRV(ierr);

	value = 1.0;
	ierr = MatSetValue(B,0,0,value,INSERT_VALUES);CHKERRV(ierr);
	ierr = MatSetValue(B,1,n-1,value,INSERT_VALUES);CHKERRV(ierr);

	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)B,"dirichlet equality");CHKERRV(ierr);
*/
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
PermonSolver<VectorBase>::~PermonSolver(){
	LOG_FUNC_BEGIN

	/* free temp vectors */
	free_temp_vectors();
	
	LOG_FUNC_END
}

/* prepare temp_vectors */
template<class VectorBase>
void PermonSolver<VectorBase>::allocate_temp_vectors(){
	LOG_FUNC_BEGIN

	GeneralVector<VectorBase> *pattern = qpdata->get_b(); /* I will allocate temp vectors subject to linear term */
	Ad = new GeneralVector<VectorBase>(*pattern);	
	
	LOG_FUNC_END
}

/* destroy temp_vectors */
template<class VectorBase>
void PermonSolver<VectorBase>::free_temp_vectors(){
	LOG_FUNC_BEGIN

	free(Ad);
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void PermonSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
/*	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	output <<  " - m:          " << m << std::endl;
	output <<  " - gamma:      " << gamma << std::endl;
	output <<  " - sigma1:     " << sigma1 << std::endl;
	output <<  " - sigma2:     " << sigma2 << std::endl;
	output <<  " - alphainit:  " << alphainit << std::endl;
*/	
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}

	output.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
/*	output_local <<  " - maxit:      " << this->maxit << std::endl;
	output_local <<  " - eps:        " << this->eps << std::endl;
	output_local <<  " - debugmode: " << this->debugmode << std::endl;

	output_local <<  " - m:          " << m << std::endl;
	output_local <<  " - gamma:      " << gamma << std::endl;
	output_local <<  " - sigma1:     " << sigma1 << std::endl;
	output_local <<  " - sigma2:     " << sigma2 << std::endl;
	output_local <<  " - alphainit:  " << alphainit << std::endl;

	output_local.synchronize();
*/
	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output_global, output_local);
		coutMaster.pop();
	}

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  " - it: " << std::setw(6) << this->it_last << ", ";
	output <<  "hess mult: " << std::setw(6) << this->hessmult_last << ", ";
	output <<  "fx: " << std::setw(10) << this->fx << ", ";	
//	output <<  "norm(gP): " << std::setw(10) << this->gP << ", ";
	output <<  "used memory: " << std::setw(6) << MemoryCheck::get_virtual() << "%" << std::endl;

	output << " - ";
//	output <<  "t_solve = " << std::setw(10) << this->timer_solve.get_value_last() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output <<  "      - fx:           " << std::setw(25) << this->fx << std::endl;
	output <<  "      - fx_control:   " << std::setw(25) << fx_quadratic+fx_linear << ", log: " << std::setw(25) << log(fx_quadratic+fx_linear)/log(10) << std::endl;
	output <<  "      - fx_linear:    " << std::setw(25) << fx_linear << ", log: " << std::setw(25) << log(fx_linear)/log(10) << std::endl;
	output <<  "      - fx_quadratic: " << std::setw(25) << fx_quadratic << ", log: " << std::setw(25) << log(fx_quadratic)/log(10) << std::endl;
//	output <<  "      - norm(gP):     " << std::setw(25) << this->gP << ", log: " << std::setw(25) << log(this->gP)/log(10) << std::endl;
	output << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - it all =        " << this->it_sum << std::endl;
	output <<  " - hessmult all =  " << this->hessmult_sum << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve =      " << this->timer_solve.get_value_sum() << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	double fx_linear, fx_quadratic;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to qpdata */
	pMatrix A = *(qpdata->get_A());
	pVector b = *(qpdata->get_b());

	/* pointer to solution */
	pVector x = *(qpdata->get_x());

	/* auxiliary vectors */
	pVector Ad = *(this->Ad); /* A*p */

	Ad = A*x;
	fx_quadratic = 0.5*dot(Ad,x);
	fx_linear = -dot(b,x);
	std::streamsize ss = std::cout.precision();

	values << std::setprecision(17);

	header << "PERMON it, ";
	values << this->it_last << ", ";

	header << "PERMON hessmult, ";
	values << this->hessmult_last << ", ";

	header << "PERMON t all, ";
	values << this->timer_solve.get_value_last() << ", ";

	header << "PERMON fx, ";
	values << this->fx << ", ";

	header << "PERMON fx_linear, ";
	values << fx_linear << ", ";

	header << "PERMON fx_quadratic, ";
	values << fx_quadratic << ", ";

	values << std::setprecision(ss);

	LOG_FUNC_END
}

template<class VectorBase>
void PermonSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	header << "PERMON_sum it, ";
	values << this->it_sum << ", ";

	header << "PERMON_sum hessmult, ";
	values << this->hessmult_sum << ", ";

	header << "PERMON_sum t all, ";
	values << this->timer_solve.get_value_sum() << ", ";

	LOG_FUNC_END
}

template<class VectorBase>
std::string PermonSolver<VectorBase>::get_name() const {
	return "PERMON";
}

/* solve the problem */
template<class VectorBase>
void PermonSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN
	
	this->timer_solve.start(); /* stop this timer in the end of solution */

	int it = 0;
	int hessmult = 0;
	double fx = std::numeric_limits<double>::max();


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


}
} /* end namespace */

#endif
