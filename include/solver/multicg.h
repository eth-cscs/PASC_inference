#ifndef PASC_MULTICGSOLVER_H
#define	PASC_MULTICGSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "solver/cgqpsolver.h"
#include "data/qpdata.h"

#include "matrix/blockdiag.h"
#include "matrix/localdense.h"

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
template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(){
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)CONSTRUCTOR" << std::endl;

	qpdata = NULL;
	
	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

}

template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(const QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;

	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

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

template<class VectorBase>
void MultiCGSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(MultiCGSolver)FUNCTION: printstatus" << std::endl;

	output <<  this->get_name() << std::endl;
	output <<  " - max(it):       " << this->it_last << std::endl;
	output <<  " - max(hessmult): " << this->hessmult_last << std::endl;
	output <<  " - max(norm(g)):  " << this->fx << std::endl;	
	output <<  " - used memory:   " << MemoryCheck::get_virtual() << "%" << std::endl;

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
	QPData<VectorBase> data_sub; /* data of inner cg solver */

	/* get number of blocks and blocks */
	BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase>> *A = dynamic_cast<BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase>> *>(qpdata->get_A());
	LocalDenseMatrix<VectorBase> **blocks = A->get_blocks();
	LocalDenseMatrix<VectorBase> *A_sub; 
	int nmb_blocks = A->get_nmb_blocks();	
	int blocksize = A->get_blocksize();

	/* get vectors */
	typedef GeneralVector<VectorBase> (&pVector);
	pVector b  = *(qpdata->get_b());
	pVector x  = *(qpdata->get_x());
	pVector x0 = *(qpdata->get_x0());

	GeneralVector<VectorBase> b_sub(blocksize);
	GeneralVector<VectorBase> x_sub(blocksize);
	GeneralVector<VectorBase> x0_sub(blocksize);

	/* through blocks */
	for(int i=0;i<nmb_blocks;i++){
		/* get data for subproblem */
		A_sub = blocks[i];
		b_sub = 1*b(i*blocksize, (i+1)*blocksize - 1); // TODO: something is wrong with petscvector
		x_sub = 1*x(i*blocksize, (i+1)*blocksize - 1); // TODO: something is wrong with petscvector
		x0_sub = 1*x0(i*blocksize, (i+1)*blocksize - 1); // TODO: something is wrong with petscvector

		/* set data of subproblem to ... subproblem data :) */
		data_sub.set_A(A_sub);
		data_sub.set_b(&b_sub);
		data_sub.set_x(&x_sub);
		data_sub.set_x0(&x0_sub);
		
		/* create new instance of solver, during the constructor the new temp vectors are initialized */
		cgsolver = new CGQPSolver<VectorBase>(data_sub);
		
		/* copy settings */
		
		/* solve this subproblem */
		cgsolver->solve();

		/* update iteration counter */
		if(cgsolver->get_fx() > this->fx) this->fx = cgsolver->get_fx();
		if(cgsolver->get_it() > this->it_last) this->it_last = cgsolver->get_it();
		if(cgsolver->get_hessmult() > this->hessmult_last) this->hessmult_last = cgsolver->get_hessmult();

		/* solution back to global vector */
		x(i*blocksize, (i+1)*blocksize - 1) = x_sub;
		
		free(cgsolver);

		// TODO: deal with iteration counters (max?)
	}


	
}

template<class VectorBase>
double MultiCGSolver<VectorBase>::get_fx() const {
	if(setting.debug_mode >= 11) coutMaster << "(MultiCGSolver)FUNCTION: get_fx()" << std::endl;
	
	return this->fx;	
}

template<class VectorBase>
int MultiCGSolver<VectorBase>::get_it() const {
	return this->it_last;
}

template<class VectorBase>
int MultiCGSolver<VectorBase>::get_hessmult() const {
	return this->hessmult_last;
}


} /* end namespace */

#endif
