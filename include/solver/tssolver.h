#ifndef PASC_TSSOLVER_H
#define	PASC_TSSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/tsdata.h"
#include "model/tsmodel.h"

#define TSSOLVER_DEFAULT_MAXIT 1000;
#define TSSOLVER_DEFAULT_EPS 0.0001;

namespace pascinference {

/* settings */
class TSSolverSetting : public GeneralSolverSetting {
	protected:
		int maxit;
		double eps;

	public:
		TSSolverSetting() {
			maxit = TSSOLVER_DEFAULT_MAXIT;
			eps = TSSOLVER_DEFAULT_EPS;
			
		};
		~TSSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << " TSSolverSettings:" << std::endl;
			output << "   - maxit: " << maxit << std::endl;
			output << "   - eps: " << eps << std::endl;

		};
		
};


/* TSSolver */ 
template<class VectorBase>
class TSSolver: public GeneralSolver {
	protected:
		const TSData<VectorBase> *data; /* data on which the solver operates */

		GeneralSolver *gammasolver; /* to solve inner gamma problem */
		GeneralSolver *thetasolver; /* to solve inner theta problem */

		TSModel<VectorBase> *model; /* pointer to used time-series model */

	public:
		TSSolverSetting setting;

		TSSolver();
		TSSolver(const TSData<VectorBase> &new_data); 
		//TODO: write constructor with given gammasolver and theta solver
		~TSSolver();

		/* with given gamma and theta solvertype */
		void solve(SolverType gammasolvertype, SolverType thetasolvertype);
		
		/* with one solvertype */
		virtual void solve(SolverType solvertype) {
			this->solve(solvertype,solvertype); 
		};
		
		/* without given solvertype */
		virtual void solve() {
			this->solve(SOLVER_AUTO); 
		};

		virtual void print(std::ostream &output) const;
		std::string get_name() const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSSolver<VectorBase>::TSSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)CONSTRUCTOR" << std::endl;
	
	data = NULL;
	model = NULL;
	gammasolver = NULL; /* in this time, we don't know how to solve the problem */
	thetasolver = NULL; /* in this time, we don't know how to solve the problem */
	
}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(const TSData<VectorBase> &new_data){
	data = &new_data;
	model = data->get_model(); 

	/* we can initialize solvers - based on model */
	model->initialize_gammasolver(&gammasolver, data);	
	model->initialize_thetasolver(&thetasolver, data);	
	
}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)DESTRUCTOR" << std::endl;

	/* destroy child solvers - based on model*/
	if(gammasolver){
		model->finalize_gammasolver(&gammasolver, data);	
	}
	if(thetasolver){
		model->finalize_thetasolver(&thetasolver, data);	
	}
	
}


/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << this->get_name() << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;
	
	/* if child solvers are specified, then print also info about it */	
	output << " Gamma Solver" << std::endl;
	if(gammasolver){
		output << *gammasolver << std::endl;
	} else {
		output << " - not set" << std::endl;
	}

	output << " Theta Solver" << std::endl;
	if(thetasolver){
		output << *thetasolver << std::endl;
	} else {
		output << " - not set" << std::endl;
	}


}

template<class VectorBase>
std::string TSSolver<VectorBase>::get_name() const {
	return "Time-Series Solver";
}

/* solve the problem */
template<class VectorBase>
void TSSolver<VectorBase>::solve(SolverType gammasolvertype, SolverType thetasolvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)FUNCTION: solve" << std::endl;

	/* the gamma or theta solver wasn't specified yet */
	if(!gammasolver || !thetasolver){
		// TODO: give error - actually, gamma and theta solvers are created during constructor with data - now we don't have data, there is nothing to solve
	}

	/* which specific solver we can use to solve the problem? */
	if(gammasolvertype == SOLVER_AUTO){
			// TODO: here write sofisticated decision tree - maybe based on MODEL?
	} 
	if(thetasolvertype == SOLVER_AUTO){
			// TODO: here write sofisticated decision tree - maybe based on MODEL?
	} 

	/* now the gammasolver and thetasolver should be specified and prepared */

	/* variables */
	double L, L_old, deltaL; /* object function value */

	/* initialize value of object function */
	L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
	
	int it; 
	
	/* main cycle */
	for(it=0;it < 1;it++){
		Message_info_value(" - it = ",it);

		/* --- COMPUTE Theta --- */
//		this->model.compute_theta(this->data.get_data_vec());
		
		/* --- COMPUTE gamma --- */
//		this->model.compute_gamma(this->data.get_data_vec());

		/* compute stopping criteria */
		L_old = L;
//		L = model.get_function_value();
		deltaL = abs(L - L_old);

		/* print info about cost function */
		Message_info_value("  - L_old       = ",L_old);
		Message_info_value("  - L           = ",L);
		Message_info_value("  - |L - L_old| = ",deltaL);

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < 0.0001){
			break;
		}
		
	}

	
}


} /* end namespace */

#endif
