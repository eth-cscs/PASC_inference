#ifndef PASC_VARXH1FEMMODEL_H
#define	PASC_VARXH1FEMMODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "model/tsmodel.h"
#include "matrix/blockdiaglaplace_explicit.h"
#include "matrix/blockdiag.h"

#include "feasibleset/simplex.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#include "vector/localvector.h"
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
	return this->dim * (1 + this->dim*this->xmem + (this->umem + 1)) * this->K;
}

/* get the appropriate length of u vector */
template<class VectorBase>
int VarxH1FEMModel<VectorBase>::get_ulength(){
	return this->umem*(this->T + this->xmem);
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
	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<VectorBase>(*gammadata->get_x0(),this->K, this->epssqr)); /* create new blockdiagonal matrix */
	gammadata->set_feasibleset(new SimplexFeasibleSet<VectorBase>(this->T,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<VectorBase>(*gammadata);
	
	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

}

/* prepare theta solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	/* in this case, theta problem is a sequence of unconstrained QP problems */

	/* prepare array with matrices */
	int nmb_blocks = this->K*this->dim;
	LocalDenseMatrix<VectorBase> **blocks = new LocalDenseMatrix<VectorBase>*[nmb_blocks];

	int blocksize = 1 + this->dim*this->xmem + (this->umem + 1); /* see get_thetavectorlength() */
	for(int i=0;i<nmb_blocks;i++){
		
	}

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
	thetadata->set_A(new BlockDiagMatrix<VectorBase,LocalDenseMatrix<VectorBase>>(nmb_blocks, blocks));


	

	/* create solver */
//	*thetasolver = new DiagSolver<VectorBase>(*thetadata);
	
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
	
	/* destroy data */
//	int k,n; /* through clusters, through dimensions */
//	for(k=0;k<this->K;k++){ // TODO: parallel
//		for(n=0;n<this->dim;n++){
//			free(thetadata[k*this->dim + n]);

			/* MATRIX */
//			if(n==0){ 
				/* if this is a first dimension, then create matrix */
//				free(thetadata[k*this->dim]->get_A());
//			}
			
			/* rhs */
//			free(thetadata[k*this->dim + n]->get_b()); 
			
			/* solution */
//			free(thetadata[k*this->dim + n]->get_x0()); // TODO: make it in a such way, that we will not need any additional storage
//			free(thetadata[k*this->dim + n]->get_x());  // TODO: make it in a such way, that we will not need any additional storage
			
//		}
		
//	}

	/* destroy solver */
//	free(*thetasolver);

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

}

/* update theta solver */
template<class VectorBase>
void VarxH1FEMModel<VectorBase>::update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	/* update theta solver - prepare new matrix vector and right-hand side vector */	

}

template<class VectorBase>
QPData<VectorBase>* VarxH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
QPData<VectorBase>* VarxH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}



} /* end namespace */

#endif
