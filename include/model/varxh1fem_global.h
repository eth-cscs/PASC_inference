#ifndef PASC_VARXH1FEMMODEL_GLOBAL_H
#define	PASC_VARXH1FEMMODEL_GLOBAL_H

#ifndef USE_PETSCVECTOR
 #error 'VARXH1FEMMODEL_GLOBAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "model/tsmodel_global.h"

/* gamma problem */
//#include "matrix/blockdiaglaplace_explicit.h"
#include "matrix/blockdiaglaplace_vector.h"

#include "feasibleset/simplex_local.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

/* theta problem */
#include "solver/multicg_global.h"
#include "matrix/blockdiag.h"

//#include "vector/localvector.h"
#include "matrix/localdense.h"

#include "data/tsdata_global.h"


namespace pascinference {



/* KMEANSH1MODEL_GLOBAL */ 
class VarxH1FEMModel_Global: public TSModel_Global {
	protected:
		QPData<PetscVector> *gammadata; /**< QP with simplex */
	 	QPData<PetscVector> *thetadata; /**< QP with blockdiag with dim*K blocks, will be solved by multiCG  */

		/* model specific variables */
		double epssqr; /**< penalty coeficient */
		int xmem; /**< size of memory for x */
		int umem; /**< size of memory for u */

	public:
		VarxH1FEMModel_Global(int T, int xdim, int K_num, int *K, int xmem, int umem, double epssqr);
		~VarxH1FEMModel_Global();

		void print(std::ostream &output) const;
		void print(std::ostream &output_global, std::ostream &output_local) const;
		std::string get_name() const;

		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData_Global *tsdata);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData_Global *tsdata);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData_Global *tsdata);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData_Global *tsdata);

		void update_gammasolver(GeneralSolver *gamma_solver, const TSData_Global *tsdata);
		void update_thetasolver(GeneralSolver *theta_solver, const TSData_Global *tsdata);
	
		void generate_data(TSData_Global *tsdata);

		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData_Global *tsdata);

		QPData<PetscVector> *get_gammadata() const;
		QPData<PetscVector> *get_thetadata() const;
		
};

} // end of namespace


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
VarxH1FEMModel_Global::VarxH1FEMModel_Global(int new_T, int new_xdim, int new_Knum, int *new_K, int new_xmem, int new_umem, double new_epssqr) {
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = new_T;
	this->Knum = new_Knum;
	this->K = new_K;

	this->Klocal = new_K[GlobalManager.get_rank()]; // TODO: Klocal will be array


	this->xdim = new_xdim;

	this->xmem = new_xmem;
	this->umem = new_umem;
	this->epssqr = new_epssqr;

	/* compute vector lengths */
	/* for data prepare vector of length T and spit it into processors */
	Vec layout;

	/* try to make a global vector of length T and then get size of local portion */
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,this->T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&(this->Tlocal)) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */
	
	this->datavectorlength_global = xdim*this->T;
	this->datavectorlength_local = xdim*this->Tlocal;

/*
	int udim;
	if(this->umem > 0){
		udim = 1;
	} else {
		udim = 0;
	}	
	this->ulength_global = udim*(this->T + this->xmem);
*/
	this->ulength_global = 0; // TODO
	this->ulength_local = 0;


	/* prepage global vectors with variables - theta and gamma */
	int Ksum = sum_array(this->Knum, this->K);
	
	this->gammavectorlength_global = Ksum*this->T;
	this->gammavectorlength_local = Klocal*this->T;
	this->thetavectorlength_global = xdim * (1 + this->xdim*this->xmem + this->umem) * Ksum; /* xdim * (mu + A + B) * K */
	this->thetavectorlength_local = xdim * (1 + this->xdim*this->xmem + this->umem) * Klocal; /* xdim * (mu + A + B) * K */
	
}

/* destructor */
VarxH1FEMModel_Global::~VarxH1FEMModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel_Global)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about model */
void VarxH1FEMModel_Global::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << this->T << std::endl;

	output <<  " - K:       ";
	print_array(output, this->Knum, this->K);
	output << std::endl;

	output <<  " - xdim:    " << this->xdim << std::endl;
	output <<  " - xmem:    " << this->xmem << std::endl;
	output <<  " - umem:    " << this->umem << std::endl;
	output <<  " - epssqr:  " << this->epssqr << std::endl;

	output <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output <<  " - ulength:     " << this->ulength_global << std::endl;
	output <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;
		
}

/* print info about model */
void VarxH1FEMModel_Global::print(std::ostream &output_global, std::ostream &output_local) const {
	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info:  " << std::endl;
	output_global <<  "  - T:       " << this->T << std::endl;

	output_global <<  "  - K:       ";
	print_array(output_global, this->Knum, this->K);
	output_global << std::endl;

	output_global <<  "  - xdim:    " << this->xdim << std::endl;
	output_global <<  "  - xmem:    " << this->xmem << std::endl;
	output_global <<  "  - umem:    " << this->umem << std::endl;
	output_global <<  "  - epssqr:  " << this->epssqr << std::endl;

	output_global <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output_global <<  " - ulength:     " << this->ulength_global << std::endl;
	output_global <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output_global <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local info:  " << std::endl;
	output_local <<   "  - T=" << this->Tlocal << ",\t K=" << this->Klocal << ",\t datalength=" << this->datavectorlength_local << ", \t ulength=" << this->ulength_local << ",\t gammalength=" << this->gammavectorlength_local << ",\t thetalength=" << this->thetavectorlength_local << std::endl;
		
}

/* get name of the model */
std::string VarxH1FEMModel_Global::get_name() const {
	return "VARX-H1-FEM Global Time-Series Model";	
}

/* prepare gamma solver */
void VarxH1FEMModel_Global::initialize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<PetscVector>();
	
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */

	gammadata->set_A(new BlockDiagLaplaceVectorMatrix<PetscVector>(*gammadata->get_x0(),this->Klocal, this->epssqr)); /* create new blockdiagonal matrix */
//	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<PetscVector>(*gammadata->get_x0(),this->K, this->epssqr)); /* create new blockdiagonal matrix */
	gammadata->set_feasibleset(new SimplexFeasibleSet_Local(this->T,this->Klocal)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<PetscVector>(*gammadata);

	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	/* set default gammasolver (qpsolver) type */
	// TODO: deprecated, decision tree implemented in QPSolver
	//dynamic_cast<QPSolver<PetscVector> *>(*gammasolver)->setting.child_solvertype = SOLVER_SPGQP;
}

/* prepare theta solver */
void VarxH1FEMModel_Global::initialize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
	/* in this case, theta problem is a sequence of unconstrained QP problems */

	/* prepare array with matrices */
	int nmb_blocks = this->Klocal*this->xdim;
	LocalDenseMatrix<PetscVector> **blocks = new LocalDenseMatrix<PetscVector>*[nmb_blocks];

	int blocksize = 1 + this->xdim*this->xmem + this->umem; /* see get_thetavectorlength() */

	int k,n; /* through clusters, through dimensions */
	for(k=0;k<this->Klocal;k++){ // TODO: parallel
		for(n=0;n<this->xdim;n++){
			if(n==0){ 
				/* if this is a first dimension, then create matrix */
				blocks[k*this->xdim + n] = new LocalDenseMatrix<PetscVector>(blocksize, blocksize);
			} else {
				/* else get matrix from first dimension */
				blocks[k*this->xdim + n] = blocks[k*this->xdim];
			}
		}
	}
	
	/* create data */
	thetadata = new QPData<PetscVector>();
	thetadata->set_x0(tsdata->get_thetavector()); /* the initial approximation of QP problem is gammavector */
	thetadata->set_x(tsdata->get_thetavector()); /* the solution of QP problem is gamma */
	thetadata->set_b(new GeneralVector<PetscVector>(*thetadata->get_x0())); /* create new linear term of QP problem */

	thetadata->set_A(new BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> >(nmb_blocks, blocks, blocksize));

	/* create solver */
	*thetasolver = new MultiCGSolver_Global(*thetadata);
	
}

/* destroy gamma solver */
void VarxH1FEMModel_Global::finalize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
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
void VarxH1FEMModel_Global::finalize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
	/* destroy blocks */

	/* prepare pointer to child, I need to change pointer from general matrix to block diag (to be able to call get_block() ) */
	BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *A = dynamic_cast<BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *>(thetadata->get_A());
	LocalDenseMatrix<PetscVector> **blocks = A->get_blocks();

	int k; /* through clusters, through dimensions */
	for(k=0;k<this->Klocal;k++){ // TODO: parallel
		/* if this is a first dimension, then destroy matrix */
		free(blocks[k*this->xdim]);
	}
	
	free(thetadata->get_A());
	
	/* destroy other data */
	free(thetadata->get_b());
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

}

double VarxH1FEMModel_Global::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData_Global *tsdata){
	
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<PetscVector> *)gammasolver)->get_fx();
}

QPData<PetscVector>* VarxH1FEMModel_Global::get_gammadata() const {
	return gammadata;
}

QPData<PetscVector>* VarxH1FEMModel_Global::get_thetadata() const {
	return thetadata;
}


void VarxH1FEMModel_Global::update_gammasolver(GeneralSolver *gammasolver, const TSData_Global *tsdata){
	/* update gamma_solver data - prepare new linear term */

	typedef GeneralVector<PetscVector> (&pVector);

	pVector M = *(tsdata->get_thetavector()); /* parameters of the model */
	pVector X = *(tsdata->get_datavector()); /* original data */
	pVector U = *(tsdata->get_u()); /* external forces */
	pVector b = *(gammadata->get_b()); /* rhs of QP gamma problem */

	/* get constants of the problem */
	int K = this->Klocal;
	int dim = this->xdim;
	int xmem = this->xmem;
	int umem = this->umem;
	int blocksize = 1 + this->xdim*this->xmem + this->umem; /* see get_thetavectorlength() */

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
				M_idxstart = k*blocksize*this->xdim + n*blocksize;
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
void VarxH1FEMModel_Global::update_thetasolver(GeneralSolver *thetasolver, const TSData_Global *tsdata){
	/* update theta solver - prepare new BlockDiag matrix data and right-hand side vector */	

	/* pointers to data */
	Vec X = tsdata->get_datavector()->get_vector();
//	Vec U = tsdata->get_u()->get_vector();
	Vec b_global = thetadata->get_b()->get_vector();
	Vec gamma_global = gammadata->get_x()->get_vector();

	/* local variables */
//	Vec b_local;
	const double *gamma_arr;
	double *b_arr;
	
	/* get constants of the problem */
	int Klocal = this->Klocal;
	int T = this->T;
	int xdim = this->xdim;
//	int xmem = this->xmem;
//	int umem = this->umem;
//	int blocksize = 1 + this->xdim*this->xmem + this->umem; /* see get_thetavectorlength() */

	/* blocks of the matrix */
	BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *A = dynamic_cast<BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *>(thetadata->get_A());
	LocalDenseMatrix<PetscVector> **blocks = A->get_blocks();
	LocalDenseMatrix<PetscVector> *block;

	/* index set with components of x corresponding to one dimension */
	IS xn_is;
	Vec xn_global; /* global vector with values of one dimension */
	Vec xn_local;  /* scattered xn_global to each processor */
	TRY( VecCreateSeq(PETSC_COMM_SELF, T, &xn_local) );
	const double *xn_local_arr; /* for read */

	VecScatter ctx; /* for scattering xn_global to xn_local */

//	TRY( VecView(X,PETSC_VIEWER_STDOUT_WORLD) );
//	TRY( VecView(b,PETSC_VIEWER_STDOUT_WORLD) );
//	TRY( VecView(gamma,PETSC_VIEWER_STDOUT_WORLD) );
	
	/* go through clusters and fill matrices */
	int k;
	int n; // n2
	int t;
	double value;


	TRY( VecGetArrayRead(gamma_global, &gamma_arr) );
	TRY( VecGetArray(b_global, &b_arr) );

	/* matrix */
	for(k=0;k<Klocal;k++){
		block = blocks[k*xdim];
			
		/* fill matrix */
		value = sum_subarray(k*T, (k+1)*T-1, gamma_arr);
		block->set_value(0, 0, value);
			
	}	

	/* rhs */
	for(n=0;n<xdim;n++){
		/* get global xn_global */
		TRY( ISCreateStride(PETSC_COMM_SELF, T, n, xdim, &xn_is) );
		TRY( VecGetSubVector(X,xn_is,&xn_global) );
		
		/* scatter xn_global to xn_local */
		TRY( VecScatterCreateToAll(xn_global,&ctx,&xn_local) );
		TRY( VecScatterBegin(ctx, xn_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
		TRY( VecScatterEnd(ctx, xn_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
		TRY( VecScatterDestroy(&ctx) );

		/* now I have my own seq vector xn_local with x(n,:) */
		TRY( VecGetArrayRead(xn_local, &xn_local_arr) );
		
		for(k=0;k<Klocal;k++){
			/* RHS - compute dot product dot( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)), gamma(k*T,(k+1)*T-1) ); */
			value = 0;
			for(t=0;t<T;t++){
				value += xn_local_arr[t]*gamma_arr[t];
			}
			b_arr[k*xdim + n] = value;
		}	

		TRY( VecGetArrayRead(xn_local, &xn_local_arr) );

		/* restore subvector with xn_global */
		TRY( VecRestoreSubVector(X,xn_is,&xn_global) );
		TRY( ISDestroy(&xn_is) );

		
	}
	TRY( VecRestoreArrayRead(gamma_global, &gamma_arr) );
	TRY( VecRestoreArray(b_global, &b_arr) );

	TRY( VecDestroy(&xn_local) );

	
//	for(k=0;k<Klocal;k++){
		/* get actual block */
//		block = blocks[k*xdim];

		/* get local gamma corresponding to this cluster */


		/* const row */
//		for(i=0; i < 1; i++){
			/* const column */
//			for(j=0; j < 1; j++){
//				value = sum_subarray(k*T, (k+1)*T-1, gamma_arr);
//				block->set_value(i, j, value);
//			}

			/* x columns */
//			for(j=1; j < xmem; j++){
				
//				for(n=0;n<dim;n++){
//					value = dot(X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j), gamma(k*T,(k+1)*T-1));
//					block->set_value(i, 1+j*dim+n, value);
//					block->set_value(1+j*dim+n, i, value);
//				}
//			}

			/* u columns */
//			if(umem > 0){
				/* there is u */
//				for(j=0; j < umem+1; j++){	// TODO: is this right?
//					value = dot(U(xmem-j, xmem+T-1-j), gamma(k*T,(k+1)*T-1));
//					block->set_value(i, 1+xmem*dim+j, value);
//					block->set_value(1+xmem*dim+j, i, value);
//				}
//			}

			/* rhs */
//			for(n=0;n<xdim;n++){
//				value = dot( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)), gamma(k*T,(k+1)*T-1) );
//				b_arr[k*blocksize*xdim + n*blocksize] = n;
//			}

			
//		}

		/* X rows */
//		for(i=1; i < xmem; i++){
//			for(n2=0;n2<dim;n2++){
//				temp = mul(X(xmem + n2*(xmem+T) - i, xmem+T-2  + n2*(xmem+T) - i), gamma(k*T,(k+1)*T-1));

				/* x columns */
//				for(j=1; j < xmem; j++){
//					for(n=0;n<dim;n++){
//						value = dot(temp, X(xmem + n*(xmem+T) - j, xmem+T-1  + n*(xmem+T) - j));
//						block->set_value(1+i*dim+n2, 1+j*dim+n, value);
//						block->set_value(1+j*dim+n, 1+i*dim+n2, value);
//					}
//				}

				/* u columns */
//				if(umem > 0){
//					for(j=0; j < umem+1; j++){	// TODO: is this right?
//						value = dot(temp, U(xmem-j, xmem+T-1-j));
//						block->set_value(1+i*dim+n2, 1+xmem*dim+j, value);
//						block->set_value(1+xmem*dim+j, 1+i*dim+n2, value);
//					}
//				}

				/* rhs */
//				for(n=0;n<dim;n++){
//					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
//					b(1 + k*blocksize*dim + n*blocksize + i*dim + n2) = value;
//				}

//			}
//		}

		/* U rows */
//		if(umem > 0){ // TODO: is this right?
//			for(i=0; i < umem+1; i++){
//				temp = mul( U(xmem-i, xmem+T-1-i), gamma(k*T,(k+1)*T-1));

				/* u columns */
//				for(j=0; j < umem+1; j++){
//					value = dot(temp, U(xmem-j, xmem+T-1-j));
//					block->set_value(1+xmem*dim+i, 1+xmem*dim+j, value);
//					block->set_value(1+xmem*dim+j, 1+xmem*dim+i, value);
//				}

				/* rhs */
//				for(n=0;n<dim;n++){
//					value = dot( temp, X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)) );
//					b(1+xmem*dim + k*blocksize*dim + n*blocksize + i) = value;
//				}
//			}
//		}

		
//	}




}



} /* end namespace */

#endif
