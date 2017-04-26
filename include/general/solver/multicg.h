/** @file multicg.h
 *  @brief Conjugate Gradient method for solving Quadratic Programs with BlockDiag matrix.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_MULTICGSOLVER_H
#define	PASC_MULTICGSOLVER_H

#include "pascinference.h"
#include "general/solver/qpsolver.h"
#include "general/solver/cgqpsolver.h"
#include "general/data/qpdata.h"

#include "general/algebra/matrix/blockdiag.h"
#include "general/algebra/matrix/localdense.h"

#define MULTICGSOLVER_DEFAULT_MAXIT 1000
#define MULTICGSOLVER_DEFAULT_EPS 0.0001
#define MULTICGSOLVER_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class MultiCGSolver
 *  \brief Conjugate Gradient method for solving Quadratic Programs with BlockDiag matrix
 *
 *  Solve each block-system using CG.
*/
template<class VectorBase>
class MultiCGSolver: public QPSolver<VectorBase> {
	protected:
		const QPData<VectorBase> *qpdata; /* data on which the solver operates, matrix has to be blogdiag */
	
	public:
		MultiCGSolver();
		MultiCGSolver(const QPData<VectorBase> &new_qpdata); 
		~MultiCGSolver();

		void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		void printstatus(ConsoleOutput &output) const;
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(){
	LOG_FUNC_BEGIN

	qpdata = NULL;
	
	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

	consoleArg.set_option_value("multicgsolver_maxit", &this->maxit, MULTICGSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("multicgsolver_eps", &this->eps, MULTICGSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("multicgsolver_debugmode", &this->debugmode, MULTICGSOLVER_DEFAULT_DEBUGMODE);

	LOG_FUNC_END
}

template<class VectorBase>
MultiCGSolver<VectorBase>::MultiCGSolver(const QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;

	this->fx = -1; /* max(norm(g)) */
	this->it_last = 0; /* max(it_block) */
	this->hessmult_last = 0; /* max(hessmult_block) */

	consoleArg.set_option_value("multicgsolver_maxit", &this->maxit, MULTICGSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("multicgsolver_eps", &this->eps, MULTICGSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("multicgsolver_debugmode", &this->debugmode, MULTICGSOLVER_DEFAULT_DEBUGMODE);

	LOG_FUNC_END
}


/* destructor */
template<class VectorBase>
MultiCGSolver<VectorBase>::~MultiCGSolver(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void MultiCGSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void MultiCGSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* print settings */
	output_local <<  " - maxit:      " << this->maxit << std::endl;
	output_local <<  " - eps:        " << this->eps << std::endl;
	output_local <<  " - debugmode: " << this->debugmode << std::endl;
	output_local.synchronize();

	/* print data */
	if(qpdata){
		output_global << "- data:" << std::endl;
		coutMaster.push();
		qpdata->print(output_global,output_local);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void MultiCGSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - max(it):       " << this->it_last << std::endl;
	output <<  " - max(hessmult): " << this->hessmult_last << std::endl;
	output <<  " - max(norm(g)):  " << this->fx << std::endl;	
	output <<  " - used memory:   " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void MultiCGSolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
std::string MultiCGSolver<VectorBase>::get_name() const {
	std::string return_value = "MultiCGSolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}


/* solve the problem */
template<class VectorBase>
void MultiCGSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	/* for each block prepare CG solver and solve the problem */
	CGQPSolver<VectorBase> *cgsolver; /* cg solver for one block */
	QPData<VectorBase> data_sub; /* data of inner cg solver */

	/* get number of blocks and blocks */
	BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *A = dynamic_cast<BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *>(qpdata->get_A());
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
	
	/* write info to log file */
	LOG_IT(this->it_last)	
	LOG_FX(this->fx)

	LOG_FUNC_END
}

template<class VectorBase>
double MultiCGSolver<VectorBase>::get_fx() const {
	if(this->debugmode >= 11) coutMaster << "(MultiCGSolver)FUNCTION: get_fx()" << std::endl;
	
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


#ifdef USE_PETSC

typedef petscvector::PetscVector PetscVector;
	
template<>
void MultiCGSolver<PetscVector>::solve() {
	LOG_FUNC_BEGIN

	/* for each block prepare CG solver and solve the problem */
	CGQPSolver<PetscVector> *cgsolver; /* cg solver for one block */
	QPData<PetscVector> data_sub; /* data of inner cg solver */

	/* get number of blocks and blocks */
	BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *A = dynamic_cast<BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *>(qpdata->get_A());
	LocalDenseMatrix<PetscVector> **blocks = A->get_blocks();
	LocalDenseMatrix<PetscVector> *A_sub; 
	int nmb_blocks = A->get_nmb_blocks();	
	int blocksize = A->get_blocksize();

	/* get vectors (assuming seq) */
	Vec b  = qpdata->get_b()->get_vector();
	Vec x  = qpdata->get_x()->get_vector();
	Vec x0 = qpdata->get_x0()->get_vector();

	Vec b_sub;
	Vec x_sub;
	Vec x0_sub;

	IS is_sub;

	GeneralVector<PetscVector> *b_subGV;
	GeneralVector<PetscVector> *x_subGV;
	GeneralVector<PetscVector> *x0_subGV;

	/* through blocks */
	for(int i=0;i<nmb_blocks;i++){
		/* get data for subproblem */
		A_sub = blocks[i];

		/* get global xn_global */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, blocksize, i*blocksize, 1, &is_sub) );
		TRYCXX( VecGetSubVector(b,is_sub,&b_sub) );
		TRYCXX( VecGetSubVector(x,is_sub,&x_sub) );
		TRYCXX( VecGetSubVector(x0,is_sub,&x0_sub) );

		/* from petscvectors to general vectors */
		b_subGV	= new GeneralVector<PetscVector>(b_sub);
		x_subGV	= new GeneralVector<PetscVector>(x_sub);
		x0_subGV = new GeneralVector<PetscVector>(x0_sub);
	
		/* set data of subproblem to ... subproblem data :) */
		data_sub.set_A(A_sub);
		data_sub.set_b(b_subGV);
		data_sub.set_x(x_subGV);
		data_sub.set_x0(x0_subGV);

		/* create new instance of solver, during the constructor the new temp vectors are initialized */
		cgsolver = new CGQPSolver<PetscVector>(data_sub);
		
		/* copy settings */


		/* solve this subproblem */
		cgsolver->solve();

		/* update iteration counter */
		if(cgsolver->get_fx() > this->fx) this->fx = cgsolver->get_fx();
		if(cgsolver->get_it() > this->it_last) this->it_last = cgsolver->get_it();
		if(cgsolver->get_hessmult() > this->hessmult_last) this->hessmult_last = cgsolver->get_hessmult();

		/* solution back to global vector */
		/* restore subvector with xn_global */
		TRYCXX( VecRestoreSubVector(b,is_sub,&b_sub) );
		TRYCXX( VecRestoreSubVector(x,is_sub,&x_sub) );
		TRYCXX( VecRestoreSubVector(x0,is_sub,&x0_sub) );
		TRYCXX( ISDestroy(&is_sub) );
		
		/* destroy globalvectors */
		free(b_subGV);
		free(x_subGV);
		free(x0_subGV);
			
		/* destroy sub solver */	
		free(cgsolver);


		// TODO: deal with iteration counters (max?)
	}

	LOG_FUNC_END
}

#endif


}
} /* end namespace */

#endif
