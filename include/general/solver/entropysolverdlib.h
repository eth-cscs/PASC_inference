/** @file entropysolverdlib.h
 *  @brief Solver which solves problem with integrals using DLib library
 *
 *  @author Anna Marchenko & Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERDLIB_H
#define	PASC_ENTROPYSOLVERDLIB_H

#include "general/solver/generalsolver.h"
#include "general/data/entropydata.h"

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
		class ExternalContent;

		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		Timer timer_solve; /**< total solution time */
		Timer timer_compute_moments; /**< time of moment computation */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments */
		GeneralVector<VectorBase> *x_power; /**< temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< temp vector for storing power of x * gamma_k */

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

/* define blank external content for general VectorBase */
template<class VectorBase>
class EntropySolverDlib<VectorBase>::ExternalContent {
};


}
} /* end namespace */

#endif
