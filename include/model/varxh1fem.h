#ifndef PASC_VARXH1FEMMODEL_H
#define	PASC_VARXH1FEMMODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "model/tsmodel.h"

/* gamma problem */
//#include "matrix/blockdiaglaplace_explicit.h"
#include "matrix/blockdiaglaplace_vector.h"

#include "feasibleset/simplex.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

/* theta problem */
#include "solver/multicg.h"
#include "matrix/blockdiag.h"

//#include "vector/localvector.h"
#include "matrix/localdense.h"

#include "data/tsdata.h"


namespace pascinference {

/* KMEANSH1MODEL */ 
template<class VectorBase>
class VarxH1FEMModel: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex */
	 	QPData<VectorBase> *thetadata; /**< QP with blockdiag with dim*K blocks, will be solved by multiCG  */

		/* model specific variables */
		double epssqr; /**< penalty coeficient */
		int xmem; /**< size of memory for x */
		int umem; /**< size of memory for u */

	public:
		VarxH1FEMModel(int T, int dimx, int K, int xmem, int umem, double epssqr);
		~VarxH1FEMModel();

		void print(std::ostream &output) const;
		std::string get_name() const;

		int get_datavectorlength();
		int get_gammavectorlength();
		int get_thetavectorlength();
		int get_ulength();
		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);

		void update_gammasolver(GeneralSolver *gamma_solver, const TSData<VectorBase> *tsdata);
		void update_thetasolver(GeneralSolver *theta_solver, const TSData<VectorBase> *tsdata);
	
		void generate_data(TSData<VectorBase> *tsdata);

		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata);

		QPData<VectorBase> *get_gammadata() const;
		QPData<VectorBase> *get_thetadata() const;
		
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
VarxH1FEMModel<VectorBase>::VarxH1FEMModel(int new_T, int new_dimx, int new_K, int new_xmem, int new_umem, double new_epssqr) {
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = new_T;
	this->K = new_K;
	this->dim = new_dimx;

	this->xmem = new_xmem;
	this->umem = new_umem;
	this->epssqr = new_epssqr;
}

/* destructor */
template<class VectorBase>
VarxH1FEMModel<VectorBase>::~VarxH1FEMModel(){
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about problem */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << this->T << std::endl;
	output <<  " - K:       " << this->K << std::endl;
	output <<  " - dimx:    " << this->dim << std::endl;
	output <<  " - xmem:    " << this->xmem << std::endl;
	output <<  " - umem:    " << this->umem << std::endl;
	output <<  " - epssqr:  " << this->epssqr << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string VarxH1FEMModel<VectorBase>::get_name() const {
	return "VARX-H1-FEM Time-Series Model";	
}

/* get the appropriate length of datavector */
template<class VectorBase>
int VarxH1FEMModel<VectorBase>::get_datavectorlength(){
	return this->dim*(this->T + this->xmem);
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int VarxH1FEMModel<VectorBase>::get_gammavectorlength(){
	return this->K*this->T;
}

/* get the appropriate length of thetavector */
/* dimx * (mu + A + B) * K */
template<class VectorBase>
int VarxH1FEMModel<VectorBase>::get_thetavectorlength(){
	return this->dim * (1 + this->dim*this->xmem + this->umem) * this->K;
}

/* get the appropriate length of u vector */
template<class VectorBase>
int VarxH1FEMModel<VectorBase>::get_ulength(){
//	return this->umem*(this->T + this->xmem); // TODO: are you sure?
	int dimu;
	if(this->umem > 0){
		dimu = 1;
	} else {
		dimu = 0;
	}	
	return dimu*(this->T + this->xmem);
}


/* prepare gamma solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<VectorBase>();
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<VectorBase>(*gammadata->get_x0())); /* create new linear term of QP problem */

	gammadata->set_A(new BlockDiagLaplaceVectorMatrix<VectorBase>(*gammadata->get_x0(),this->K, this->epssqr)); /* create new blockdiagonal matrix */
//	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<VectorBase>(*gammadata->get_x0(),this->K, this->epssqr)); /* create new blockdiagonal matrix */

	gammadata->set_feasibleset(new SimplexFeasibleSet<VectorBase>(this->T,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<VectorBase>(*gammadata);
	
	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	/* set default gammasolver (qpsolver) type */
	// TODO: deprecated, decision tree implemented in QPSolver
	//dynamic_cast<QPSolver<VectorBase> *>(*gammasolver)->setting.child_solvertype = SOLVER_SPGQP;
}

/* prepare theta solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	/* in this case, theta problem is a sequence of unconstrained QP problems */

	/* prepare array with matrices */
	int nmb_blocks = this->K*this->dim;
	LocalDenseMatrix<VectorBase> **blocks = new LocalDenseMatrix<VectorBase>*[nmb_blocks];

	int blocksize = 1 + this->dim*this->xmem + this->umem; /* see get_thetavectorlength() */

	int k,n; /* through clusters, through dimensions */
	for(k=0;k<this->K;k++){ // TODO: parallel
		for(n=0;n<this->dim;n++){
			if(n==0){ 
				/* if this is a first dimension, then create matrix */
				blocks[k*this->dim + n] = new LocalDenseMatrix<VectorBase>(blocksize, blocksize);
			} else {
				/* else get matrix from first dimension */
				blocks[k*this->dim + n] = blocks[k*this->dim];
			}
		}
	}
	
	/* create data */
	thetadata = new QPData<VectorBase>();
	thetadata->set_x0(tsdata->get_thetavector()); /* the initial approximation of QP problem is gammavector */
	thetadata->set_x(tsdata->get_thetavector()); /* the solution of QP problem is gamma */
	thetadata->set_b(new GeneralVector<VectorBase>(*thetadata->get_x0())); /* create new linear term of QP problem */

	thetadata->set_A(new BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> >(nmb_blocks, blocks, blocksize));

	/* create solver */
	*thetasolver = new MultiCGSolver<VectorBase>(*thetadata);
	
}

/* destroy gamma solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	/* I created this objects, I should destroy them */

	/* destroy data */
	free(gammadata->get_b());
	free(gammadata->get_A());
	free(gammadata->get_feasibleset());
	free(gammadata);

	/* destroy solver */
	free(*gammasolver);
	
}

/* destroy theta solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	/* destroy blocks */

	/* prepare pointer to child, I need to change pointer from general matrix to block diag (to be able to call get_block() ) */
	BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *A = dynamic_cast<BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *>(thetadata->get_A());
	LocalDenseMatrix<VectorBase> **blocks = A->get_blocks();

	int k; /* through clusters, through dimensions */
	for(k=0;k<this->K;k++){ // TODO: parallel
		/* if this is a first dimension, then destroy matrix */
		free(blocks[k*this->dim]);
	}
	
	free(thetadata->get_A());
	
	/* destroy other data */
	free(thetadata->get_b());
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

}

template<class VectorBase>
double VarxH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

/* ---------------------- GENERAL ----------------------- */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::update_gammasolver(GeneralSolver *gammasolver, const TSData<VectorBase> *tsdata){
	/* update gamma_solver data - prepare new linear term */

	typedef GeneralVector<VectorBase> (&pVector);

	pVector M = *(tsdata->get_thetavector()); /* parameters of the model */
	pVector X = *(tsdata->get_datavector()); /* original data */
	pVector U = *(tsdata->get_u()); /* external forces */
	pVector b = *(gammadata->get_b()); /* rhs of QP gamma problem */

	/* get constants of the problem */
	int K = this->K;
	int dim = this->dim;
	int xmem = this->xmem;
	int umem = this->umem;
	int blocksize = 1 + this->dim*this->xmem + this->umem; /* see get_thetavectorlength() */

	int t,n,k,j,n2;
	int M_idxstart;
	int Xdimlength = X.size()/(double)dim;
	int T = Xdimlength - xmem;
	double dot_sum, value, valuemodel;
	
	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			dot_sum = 0.0;
			for(n=0;n<dim;n++){
				value = X(n*Xdimlength+xmem+t); 
				
				/* compute model_value */
				M_idxstart = k*blocksize*this->dim + n*blocksize;
				valuemodel = M(M_idxstart + 0); /* mu_K,n */
				for(j=1;j<xmem;j++){ /* add A_j*X(t-j) */
					for(n2=0;n2<dim;n2++){
						valuemodel += M(M_idxstart + 1 + j*dim + n2)*X(Xdimlength*n2 + xmem - j + t);
					}
				}
				
				if(umem > 0){
					/* there is u */
					for(j=0;j<umem-1;j++){ /* add B_j*U(t) */
						valuemodel += M(M_idxstart + 1 + xmem*dim + j)*U(xmem - j + t);
					}
				}
				
				value = value - valuemodel; //TODO: get in MinLin?
				dot_sum += value*value;
			}
			b(k*T + t) = -(dot_sum);
		}
	}	

}

/* update theta solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	/* update theta solver - prepare new BlockDiag matrix data and right-hand side vector */	

	/* pointers to data */
	typedef GeneralVector<VectorBase> (&pVector);
	pVector X = *(tsdata->get_datavector());
	pVector U = *(tsdata->get_u());
	pVector b = *(thetadata->get_b());
	pVector gamma = *(gammadata->get_x());

	/* get constants of the problem */
	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	int xmem = this->xmem;
	int umem = this->umem;
	int blocksize = 1 + this->dim*this->xmem + this->umem; /* see get_thetavectorlength() */

	/* auxiliary vector */
	VectorBase temp(T);

	/* blocks of the matrix */
	BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *A = dynamic_cast<BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase> > *>(thetadata->get_A());
	LocalDenseMatrix<VectorBase> **blocks = A->get_blocks();
	LocalDenseMatrix<VectorBase> *block;
	
	/* go through clusters and fill matrices */
	int k,n,n2;
	int i,j;
	double value;
	for(k=0;k<K;k++){
		/* get actual block */
		block = blocks[k*dim];
		
		/* const row */
		for(i=0; i < 1; i++){
			/* const column */
			for(j=0; j < 1; j++){
//				value = T;
				value = sum(gamma(k*T,(k+1)*T-1));
				block->set_value(i, j, value);
			}

			/* x columns */
			for(j=1; j < xmem; j++){
				
				for(n=0;n<dim;n++){
//					value = sum(X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
//					value = dot(X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j), gamma(k*T,(k+1)*T-1));
					value = dot(X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j), gamma(k*T,(k+1)*T-1));
					block->set_value(i, 1+j*dim+n, value);
					block->set_value(1+j*dim+n, i, value);
				}
			}

			/* u columns */
			if(umem > 0){
				/* there is u */
				for(j=0; j < umem+1; j++){	// TODO: is this right?
//					value = sum(U(xmem-j, xmem+T-1-j));
//					value = dot(U(xmem-j, xmem+T-1-j), gamma(k*T,(k+1)*T-1));
					value = dot(U(xmem-j, xmem+T-1-j), gamma(k*T,(k+1)*T-1));
					block->set_value(i, 1+xmem*dim+j, value);
					block->set_value(1+xmem*dim+j, i, value);
				}
			}

			/* rhs */
			for(n=0;n<dim;n++){
//				value = sum( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
				value = dot( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)), gamma(k*T,(k+1)*T-1) );
				b(k*blocksize*dim + n*blocksize) = value;
			}

			
		}

		/* X rows */
		for(i=1; i < xmem; i++){
			for(n2=0;n2<dim;n2++){
//				temp = mul(X(xmem-1 + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), gamma(k*T,(k+1)*T-1));
				temp = mul(X(xmem + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), gamma(k*T,(k+1)*T-1));

				/* x columns */
				for(j=1; j < xmem; j++){
					for(n=0;n<dim;n++){
//						value = dot(X(xmem-1 + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
//						value = dot(temp, X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
						value = dot(temp, X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j));
						block->set_value(1+i*dim+n2, 1+j*dim+n, value);
						block->set_value(1+j*dim+n, 1+i*dim+n2, value);
					}
				}

				/* u columns */
				if(umem > 0){
					for(j=0; j < umem+1; j++){	// TODO: is this right?
						value = dot(temp, U(xmem-j, xmem+T-1-j));
						block->set_value(1+i*dim+n2, 1+xmem*dim+j, value);
						block->set_value(1+xmem*dim+j, 1+i*dim+n2, value);
					}
				}

				/* rhs */
				for(n=0;n<dim;n++){
					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					b(1 + k*blocksize*dim + n*blocksize + i*dim + n2) = value;
				}

			}
		}

		/* U rows */
		if(umem > 0){ // TODO: is this right?
			for(i=0; i < umem+1; i++){
				temp = mul( U(xmem-i, xmem+T-1-i), gamma(k*T,(k+1)*T-1));

				/* u columns */
				for(j=0; j < umem+1; j++){
//					value = dot(U(xmem-i, xmem+T-1-i), U(xmem-j, xmem+T-1-j));				
					value = dot(temp, U(xmem-j, xmem+T-1-j));
					block->set_value(1+xmem*dim+i, 1+xmem*dim+j, value);
					block->set_value(1+xmem*dim+j, 1+xmem*dim+i, value);
				}

				/* rhs */
				for(n=0;n<dim;n++){
//					value = dot( U(xmem-i, xmem+T-1-i), X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					b(1+xmem*dim + k*blocksize*dim + n*blocksize + i) = value;
				}
			}
		}

		
	}

}

template<class VectorBase>
QPData<VectorBase>* VarxH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
QPData<VectorBase>* VarxH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}


/* PETSC SPECIFIC */
#ifdef USE_PETSCVECTOR

typedef petscvector::PetscVector GlobalPetscVector;

/* ---------------------- GENERAL ----------------------- */
template<>
void VarxH1FEMModel<GlobalPetscVector>::update_gammasolver(GeneralSolver *gammasolver, const TSData<GlobalPetscVector> *tsdata){
	/* update gamma_solver data - prepare new linear term */

	typedef GeneralVector<GlobalPetscVector> (&pVector);

	pVector M = *(tsdata->get_thetavector()); /* parameters of the model */
	pVector X = *(tsdata->get_datavector()); /* original data */
	pVector U = *(tsdata->get_u()); /* external forces */
	pVector b = *(gammadata->get_b()); /* rhs of QP gamma problem */

	/* get constants of the problem */
	int K = this->K;
	int dim = this->dim;
	int xmem = this->xmem;
	int umem = this->umem;
	int blocksize = 1 + this->dim*this->xmem + this->umem; /* see get_thetavectorlength() */

	int t,n,k,j,n2;
	int M_idxstart;
	int Xdimlength = X.size()/(double)dim;
	int T = Xdimlength - xmem;
	double dot_sum, value, valuemodel;
	
	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			dot_sum = 0.0;
			for(n=0;n<dim;n++){
				value = X.get(n*Xdimlength+xmem+t); 
				
				/* compute model_value */
				M_idxstart = k*blocksize*this->dim + n*blocksize;
				valuemodel = M.get(M_idxstart + 0); /* mu_K,n */
				for(j=1;j<xmem;j++){ /* add A_j*X(t-j) */
					for(n2=0;n2<dim;n2++){
						valuemodel += M.get(M_idxstart + 1 + j*dim + n2)*X.get(Xdimlength*n2 + xmem - j + t);
					}
				}
				
				if(umem > 0){
					/* there is u */
					for(j=0;j<umem-1;j++){ /* add B_j*U(t) */
						valuemodel += M.get(M_idxstart + 1 + xmem*dim + j)*U.get(xmem - j + t);
					}
				}
				
				value = value - valuemodel; //TODO: get in MinLin?
				dot_sum += value*value;
			}
			b(k*T + t) = -(dot_sum);
		}
	}	

}

/* update theta solver */
template<>
void VarxH1FEMModel<GlobalPetscVector>::update_thetasolver(GeneralSolver *thetasolver, const TSData<GlobalPetscVector> *tsdata){
	/* update theta solver - prepare new BlockDiag matrix data and right-hand side vector */	

	/* pointers to data */
	typedef GeneralVector<GlobalPetscVector> (&pVector);
	pVector X = *(tsdata->get_datavector());
	pVector U = *(tsdata->get_u());
	pVector b = *(thetadata->get_b());
	pVector gamma = *(gammadata->get_x());

	/* get constants of the problem */
	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	int xmem = this->xmem;
	int umem = this->umem;
	int blocksize = 1 + this->dim*this->xmem + this->umem; /* see get_thetavectorlength() */

	/* auxiliary vector */
	GlobalPetscVector temp(T);

	/* blocks of the matrix */
	BlockDiagMatrix<GlobalPetscVector,LocalDenseMatrix<GlobalPetscVector> > *A = dynamic_cast<BlockDiagMatrix<GlobalPetscVector,LocalDenseMatrix<GlobalPetscVector> > *>(thetadata->get_A());
	LocalDenseMatrix<GlobalPetscVector> **blocks = A->get_blocks();
	LocalDenseMatrix<GlobalPetscVector> *block;
	
	/* go through clusters and fill matrices */
	int k,n,n2;
	int i,j;
	double value;
	for(k=0;k<K;k++){
		/* get actual block */
		block = blocks[k*dim];
		
		/* const row */
		for(i=0; i < 1; i++){
			/* const column */
			for(j=0; j < 1; j++){
//				value = T;
				value = sum(gamma(k*T,(k+1)*T-1));
				block->set_value(i, j, value);
			}

			/* x columns */
			for(j=1; j < xmem; j++){
				
				for(n=0;n<dim;n++){
//					value = sum(X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
//					value = dot(X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j), gamma(k*T,(k+1)*T-1));
					value = dot(X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j), gamma(k*T,(k+1)*T-1));
					block->set_value(i, 1+j*dim+n, value);
					block->set_value(1+j*dim+n, i, value);
				}
			}

			/* u columns */
			if(umem > 0){
				/* there is u */
				for(j=0; j < umem+1; j++){	// TODO: is this right?
//					value = sum(U(xmem-j, xmem+T-1-j));
//					value = dot(U(xmem-j, xmem+T-1-j), gamma(k*T,(k+1)*T-1));
					value = dot(U(xmem-j, xmem+T-1-j), gamma(k*T,(k+1)*T-1));
					block->set_value(i, 1+xmem*dim+j, value);
					block->set_value(1+xmem*dim+j, i, value);
				}
			}

			/* rhs */
			for(n=0;n<dim;n++){
//				value = sum( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
				value = dot( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)), gamma(k*T,(k+1)*T-1) );
				b(k*blocksize*dim + n*blocksize) = value;
			}

			
		}

		/* X rows */
		for(i=1; i < xmem; i++){
			for(n2=0;n2<dim;n2++){
//				temp = mul(X(xmem-1 + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), gamma(k*T,(k+1)*T-1));
				temp = mul(X(xmem + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), gamma(k*T,(k+1)*T-1));

				/* x columns */
				for(j=1; j < xmem; j++){
					for(n=0;n<dim;n++){
//						value = dot(X(xmem-1 + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
//						value = dot(temp, X(xmem-1 + n*(xmem+T) - j, xmem+T-2  + n*(xmem+T) - j));
						value = dot(temp, X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j));
						block->set_value(1+i*dim+n2, 1+j*dim+n, value);
						block->set_value(1+j*dim+n, 1+i*dim+n2, value);
					}
				}

				/* u columns */
				if(umem > 0){
					for(j=0; j < umem+1; j++){	// TODO: is this right?
						value = dot(temp, U(xmem-j, xmem+T-1-j));
						block->set_value(1+i*dim+n2, 1+xmem*dim+j, value);
						block->set_value(1+xmem*dim+j, 1+i*dim+n2, value);
					}
				}

				/* rhs */
				for(n=0;n<dim;n++){
					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					b(1 + k*blocksize*dim + n*blocksize + i*dim + n2) = value;
				}

			}
		}

		/* U rows */
		if(umem > 0){ // TODO: is this right?
			for(i=0; i < umem+1; i++){
				temp = mul( U(xmem-i, xmem+T-1-i), gamma(k*T,(k+1)*T-1));

				/* u columns */
				for(j=0; j < umem+1; j++){
//					value = dot(U(xmem-i, xmem+T-1-i), U(xmem-j, xmem+T-1-j));				
					value = dot(temp, U(xmem-j, xmem+T-1-j));
					block->set_value(1+xmem*dim+i, 1+xmem*dim+j, value);
					block->set_value(1+xmem*dim+j, 1+xmem*dim+i, value);
				}

				/* rhs */
				for(n=0;n<dim;n++){
//					value = dot( U(xmem-i, xmem+T-1-i), X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
					b(1+xmem*dim + k*blocksize*dim + n*blocksize + i) = value;
				}
			}
		}

		
	}

}

#endif



} /* end namespace */

#endif
