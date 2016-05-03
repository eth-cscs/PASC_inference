#ifndef PASC_MULTICGSOLVER_GLOBAL_H
#define	PASC_MULTICGSOLVER_GLOBAL_H

#ifndef USE_PETSCVECTOR
 #error 'MULTICGSOLVER_GLOBAL_GLOBAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>

#include "solver/multicg.h"
#include "solver/qpsolver.h"
#include "solver/cgqpsolver.h"
#include "data/qpdata.h"

#include "matrix/blockdiag.h"
#include "matrix/localdense.h"

#define MULTICGSOLVER_GLOBAL_DEFAULT_MAXIT 1000;
#define MULTICGSOLVER_GLOBAL_DEFAULT_EPS 0.0001;
#define MULTICGSOLVER_GLOBAL_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class MultiCGSolver_GlobalSetting_Global : public QPSolverSetting {
	public:
		MultiCGSolver_GlobalSetting_Global() {
			this->maxit = MULTICGSOLVER_GLOBAL_DEFAULT_MAXIT;
			this->eps = MULTICGSOLVER_GLOBAL_DEFAULT_EPS;
			this->debug_mode = MULTICGSOLVER_GLOBAL_DEFAULT_DEBUG_MODE;
		};
		~MultiCGSolver_GlobalSetting_Global() {};

		virtual void print(std::ostream &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - maxit:      " << this->maxit << std::endl;
			output <<  " - eps:        " << this->eps << std::endl;
			output <<  " - debug_mode: " << this->debug_mode << std::endl;

		};

		std::string get_name() const {
			return "MultiCG_Global SolverSetting";
		};
		
};


/* MultiCGSolver_Global */ 
class MultiCGSolver_Global: public QPSolver<PetscVector> {
	protected:
		const QPData<PetscVector> *qpdata; /* data on which the solver operates, matrix has to be blogdiag */
	
	public:
		MultiCGSolver_GlobalSetting_Global setting;

		MultiCGSolver_Global();
		MultiCGSolver_Global(const QPData<PetscVector> &new_qpdata); 
		~MultiCGSolver_Global();

		void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		void print(std::ostream &output) const;
		void printstatus(std::ostream &output) const;
		void printcontent(std::ostream &output) const;
		std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
MultiCGSolver_Global::MultiCGSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

}

MultiCGSolver_Global::MultiCGSolver_Global(const QPData<PetscVector> &new_qpdata){
	qpdata = &new_qpdata;

	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

}


/* destructor */
MultiCGSolver_Global::~MultiCGSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)DESTRUCTOR" << std::endl;

}

/* print info about problem */
void MultiCGSolver_Global::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)FUNCTION: print" << std::endl;

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

void MultiCGSolver_Global::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)FUNCTION: printstatus" << std::endl;

	output <<  this->get_name() << std::endl;
	output <<  " - max(it):       " << this->it_last << std::endl;
	output <<  " - max(hessmult): " << this->hessmult_last << std::endl;
	output <<  " - max(norm(g)):  " << this->fx << std::endl;	
	output <<  " - used memory:   " << MemoryCheck::get_virtual() << "%" << std::endl;

}

/* print content of solver */
void MultiCGSolver_Global::printcontent(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)FUNCTION: printcontent" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
}

std::string MultiCGSolver_Global::get_name() const {
	return "MultiCG_Global method for QP with BlockDiag system matrix";
}

double MultiCGSolver_Global::get_fx() const {
	if(setting.debug_mode >= 11) coutMaster << "(MultiCGSolver_Global)FUNCTION: get_fx()" << std::endl;
	
	return this->fx;	
}

int MultiCGSolver_Global::get_it() const {
	return this->it_last;
}

int MultiCGSolver_Global::get_hessmult() const {
	return this->hessmult_last;
}


/* solve with Petscvector on every processor */
void MultiCGSolver_Global::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver_Global)FUNCTION: solve" << std::endl;

	/* for each block prepare CG solver and solve the problem */
	MultiCGSolver<PetscVector> *solver_local; /* local solver */

	/* global data */
	Vec x_global = qpdata->get_x()->get_vector();
	Vec x0_global = qpdata->get_x0()->get_vector();
	Vec b_global = qpdata->get_b()->get_vector();

	/* local data */
	int local_size;
	Vec x_local;
	Vec x0_local;
	Vec b_local;

	/* allocate local data */
	TRY( VecGetLocalSize(x_global, &local_size) );
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x0_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &b_local) );	

	/* get local vectors */
	TRY( VecGetLocalVector(x_global,x_local) );
	TRY( VecGetLocalVectorRead(x0_global,x0_local) );
	TRY( VecGetLocalVectorRead(b_global,b_local) );

	GeneralVector<PetscVector> x(x_local);
	GeneralVector<PetscVector> x0(x0_local);
	GeneralVector<PetscVector> b(b_local);

	/* get local vectors and prepare local data */
	QPData<PetscVector> data_local; /* data of inner cg solver */
	data_local.set_A(qpdata->get_A());
	data_local.set_b(&b);
	data_local.set_x(&x);
	data_local.set_x0(&x0);

	/* create new instance of local solver */
	solver_local = new MultiCGSolver<PetscVector>(data_local);

	/* copy settings */
		
	/* solve local problem */
	solver_local->solve();

	/* restore global vectors */
	TRY( VecRestoreLocalVector(x_global,x_local) );
	TRY( VecRestoreLocalVectorRead(x0_global,x0_local) );
	TRY( VecRestoreLocalVectorRead(b_global,b_local) );

	/* destroy local storage */
//	TRY( VecDestroy(&x_local) );	
//	TRY( VecDestroy(&x0_local) );	
//	TRY( VecDestroy(&b_local) );	

	/* free local solver */
	free(solver_local); 
	
	//TODO:temp
//	TRY( VecView(x_global,PETSC_VIEWER_STDOUT_WORLD) );
}



} /* end namespace */

#endif
