/** @file simplesolver.h
 *  @brief Solver which doesn't do anything
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_SIMPLESOLVER_H
#define	PASC_SIMPLESOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUGMODE;

#include "pascinference.h"
#include "data/simpledata.h"

namespace pascinference {
namespace solver {

/** \class SimpleSolver
 *  \brief Solver which doesn't do anything.
 *
 *  Operates on SimpleData.
*/
template<class VectorBase>
class SimpleSolver: public GeneralSolver {
	protected:
		Timer timer_solve; /**< total solution time of Simple "algorithm" */

		SimpleData<VectorBase> *simpledata; /**< data on which the solver operates */

	public:

		SimpleSolver();
		SimpleSolver(SimpleData<VectorBase> &new_simpledata); 
		~SimpleSolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printcontent(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		SimpleData<VectorBase> *get_data() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
SimpleSolver<VectorBase>::SimpleSolver(){
	LOG_FUNC_BEGIN
	
	simpledata = NULL;
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	

	LOG_FUNC_END
}

template<class VectorBase>
SimpleSolver<VectorBase>::SimpleSolver(SimpleData<VectorBase> &new_simpledata){
	LOG_FUNC_BEGIN

	simpledata = &new_simpledata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
SimpleSolver<VectorBase>::~SimpleSolver(){
	LOG_FUNC_BEGIN
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void SimpleSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << "\n";
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << "\n";
	output <<  " - eps:        " << this->eps << "\n";
	output <<  " - debugmode: " << this->debugmode << "\n";

	/* print data */
	if(simpledata){
		coutMaster.push();
		simpledata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void SimpleSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << "\n";
	
	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << "\n";
	output_global <<  " - eps:        " << this->eps << "\n";
	output_global <<  " - debugmode: " << this->debugmode << "\n";

	/* print data */
	if(simpledata){
		output_global << "- data:\n";
		coutMaster.push();
		simpledata->print(output_global,output_local);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void SimpleSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << "\n";
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%\n";

	LOG_FUNC_END
}

template<class VectorBase>
void SimpleSolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void SimpleSolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << "\n";
	
	/* print content of data */
	if(simpledata){
		output << "- data:\n";
		coutMaster.push();
		simpledata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void SimpleSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << "\n";
	output <<  " - timers\n";
	output <<  "  - t_solve =  " << this->timer_solve.get_value_sum() << "\n";

	LOG_FUNC_END
}

template<class VectorBase>
std::string SimpleSolver<VectorBase>::get_name() const {
	return "SimpleSolver";
}

template<class VectorBase>
SimpleData<VectorBase> *SimpleSolver<VectorBase>::get_data() const {
	return simpledata;
}

template<class VectorBase>
void SimpleSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->timer_solve.start(); 

	/* this solve does not anything */
	
	this->timer_solve.stop(); 

	LOG_FUNC_END
}


}
} /* end namespace */

#endif
