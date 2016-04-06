#ifndef PASC_MULTICGSOLVER_H
#define	PASC_MULTICGSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "solver/cgqpsolver.h"
#include "data/qpdata.h"

#include "matrix/blockdiag.h"

#define MULTICGSOLVER_DEFAULT_MAXIT 1000;
#define MULTICGSOLVER_DEFAULT_EPS 0.0001;
#define MULTICGSOLVER_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class MultiCGSolverSetting : public QPSolverSetting {
	public:
		MultiCGSolverSetting() {
			this->maxit = MULTICGSOLVER_DEFAULT_MAXIT;
			this->eps = MULTICGSOLVER_DEFAULT_EPS;
			this->debug_mode = MULTICGSOLVER_DEFAULT_DEBUG_MODE;
		};
		~MultiCGSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - maxit:      " << this->maxit << std::endl;
			output <<  " - eps:        " << this->eps << std::endl;
			output <<  " - debug_mode: " << this->debug_mode << std::endl;

		};

		std::string get_name() const {
			return "MultiCG SolverSetting";
		};
		
};


/* MultiCGSolver */ 
template<class VectorBase>
class MultiCGSolver: public QPSolver<VectorBase> {
	protected:
		const QPData<VectorBase> *qpdata; /* data on which the solver operates, matrix has to be blogdiag */
	
	public:
		MultiCGSolverSetting setting;

		MultiCGSolver();
		MultiCGSolver(const QPData<VectorBase> &new_qpdata); 
		~MultiCGSolver();

		void solve();

		void print(std::ostream &output) const;
		void printcontent(std::ostream &output) const;
		std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(){
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	this->fx = std::numeric_limits<double>::max();
}

template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(const QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;

	this->fx = std::numeric_limits<double>::max();
}


/* destructor */
template<class VectorBase>
MultiCGSolver<VectorBase>::~MultiCGSolver(){
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)DESTRUCTOR" << std::endl;

}

/* print info about problem */
template<class VectorBase>
void MultiCGSolver<VectorBase>::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)FUNCTION: print" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	coutMaster.push();
	setting.print(output);
	coutMaster.pop();

	/* print data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}
		
}

/* print content of solver */
template<class VectorBase>
void MultiCGSolver<VectorBase>::printcontent(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)FUNCTION: printcontent" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
}

template<class VectorBase>
std::string MultiCGSolver<VectorBase>::get_name() const {
	return "MultiCG method for QP with BlockDiag system matrix";
}


/* solve the problem */
template<class VectorBase>
void MultiCGSolver<VectorBase>::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)FUNCTION: solve" << std::endl;

	/* for each block prepare CG solver and solve the problem */
	CGQPSolver<VectorBase> *cgsolver; /* cg solver for one block */
	QPData<VectorBase> *cgsolver_data; /* data of inner cg solver */


	
}


} /* end namespace */

#endif
