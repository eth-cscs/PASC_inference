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
#include "solver/qpsolver_global.h"
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
		double *epssqr; /**< penalty coeficient */
		double epssqrlocal;
		
		int *xmem; /**< size of memory for x */
		int xmemlocal;

//		int ulength_global;
//		int ulength_local;

		void scatter_xn(Vec &x_global, Vec &xn_local, int xdim, int n);
		
	public:
		VarxH1FEMModel_Global(int T, int xdim, int num, int *K, int *xmem, double *epssqr);
		~VarxH1FEMModel_Global();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
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
VarxH1FEMModel_Global::VarxH1FEMModel_Global(int new_T, int new_xdim, int new_num, int *new_K, int *new_xmem, double *new_epssqr) {
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = new_T;
	this->xdim = new_xdim;

	this->num = new_num;
	this->K = new_K;
	this->xmem = new_xmem;
	this->epssqr = new_epssqr;

	int my_rank = GlobalManager.get_rank();

	this->Klocal = new_K[my_rank];
	this->xmemlocal = new_xmem[my_rank]; 
	this->epssqrlocal = new_epssqr[my_rank]; 

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

	/* prepage global vectors with variables - theta and gamma */
	int Ksum = sum_array(this->num, this->K);
	
	int Kxmem[this->num];
	mult_pw_array(this->num, Kxmem, this->K, this->xmem);
	int Kxmemsum = sum_array(this->num,Kxmem);
		
	this->gammavectorlength_global = Ksum*this->T - Kxmemsum;
	this->gammavectorlength_local = Klocal*(this->T - this->xmemlocal);

	this->thetavectorlength_global = xdim*(Ksum + xdim*Kxmemsum); /* all(mu + A)  */
	this->thetavectorlength_local = xdim * (1 + this->xdim*this->xmemlocal) * Klocal; /* xdim * (mu + A) * K */
	
}

/* destructor */
VarxH1FEMModel_Global::~VarxH1FEMModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(VarxH1FEMModel_Global)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about model */
void VarxH1FEMModel_Global::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K:       ";
	print_array(output, this->num, this->K);
	output << std::endl;

	output <<  " - xmem:    ";
	print_array(output, this->num, this->xmem);
	output << std::endl;

	output <<  " - epssqr:  ";
	print_array(output, this->num, this->epssqr);
	output << std::endl;

	output <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;
	
	output.synchronize();	
}

/* print info about model */
void VarxH1FEMModel_Global::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info:  " << std::endl;
	output_global <<  "  - T:       " << this->T << std::endl;
	output_global <<  "  - xdim:    " << this->xdim << std::endl;

	output_global <<  "  - K:       ";
	print_array(output_global, this->num, this->K);
	output_global << std::endl;

	output_global <<  "  - xmem:    ";
	print_array(output_global, this->num, this->xmem);
	output_global << std::endl;

	output_global <<  "  - epssqr:  ";
	print_array(output_global, this->num, this->epssqr);
	output_global << std::endl;


	output_global <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output_global <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output_global <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables:  " << std::endl;
	output_global.push();
	output_local << "T=" << std::setw(6) << this->Tlocal << ", ";
	output_local << "K=" << std::setw(6) << this->Klocal << ", ";
	output_local << "xmem=" << std::setw(6) << this->xmemlocal << ", ";
	output_local << "epssqr=" << std::setw(6) << this->epssqrlocal << ", ";
	output_local << "datalength=" << std::setw(6) << this->datavectorlength_local << ", ";
	output_local << "gammalength=" << std::setw(6) << this->gammavectorlength_local << ", ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

}

/* print model solution */
void VarxH1FEMModel_Global::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  "theta:" << std::endl;

	std::ostringstream temp;
	
	double *theta;
	TRY( VecGetArray(thetadata->get_x()->get_vector(),&theta) );
	
	int k,i,n,n2;
	int blocksize = 1 + this->xdim*this->xmemlocal;

	output_local.push();
	output_local << "- proc: " << GlobalManager.get_rank() << std::endl;
	output_local.push();
	for(k=0;k<this->Klocal;k++){
		output_local <<  "- k = " << k << std::endl;

		/* mu */
		output_local.push();
		output_local <<  "- mu = [";
		for(n=0;n<xdim;n++){
			temp << theta[k*xdim*blocksize + n*blocksize];
			output_local << temp.str();
			if(n < xdim-1){
				output_local << ", ";
			}
			temp.str("");
		}
		output_local <<  "]" << std::endl;
		output_local.pop();

		/* A */
		output_local.push();
		for(i=0;i<this->xmemlocal;i++){
			output_local <<  "- A" << i << " = [" <<  std::endl;
			output_local.push();
			for(n=0;n<xdim;n++){
				for(n2=0;n2<xdim;n2++){
					temp << theta[k*xdim*blocksize + i*xdim + n*blocksize + 1 + n2];
					output_local << temp.str() << ", ";
					temp.str("");
				}
				output_local << std::endl;
			}
			output_local.pop();
			output_local <<  "  ]" << std::endl;
		}
		output_local.pop();

	}
	output_local.pop();
	output_local.pop();

	output_local.synchronize();
	output_global.synchronize();

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

//	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<PetscVector>(*gammadata->get_x0(),this->Klocal, this->epssqr)); /* create new blockdiagonal matrix */
	gammadata->set_A(new BlockDiagLaplaceVectorMatrix<PetscVector>(*gammadata->get_x0(),this->Klocal, this->T - this->xmemlocal,this->epssqrlocal)); /* create new blockdiagonal matrix */
//	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<PetscVector>(*gammadata->get_x0(),this->K, this->epssqr)); /* create new blockdiagonal matrix */
	gammadata->set_feasibleset(new SimplexFeasibleSet_Local(this->T - this->xmemlocal,this->Klocal)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver_Global(*gammadata);

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

	int blocksize = 1 + this->xdim*this->xmemlocal; /* see get_thetavectorlength() */

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
	return ((QPSolver_Global *)gammasolver)->get_fx();
}

QPData<PetscVector>* VarxH1FEMModel_Global::get_gammadata() const {
	return gammadata;
}

QPData<PetscVector>* VarxH1FEMModel_Global::get_thetadata() const {
	return thetadata;
}

void VarxH1FEMModel_Global::update_gammasolver(GeneralSolver *gammasolver, const TSData_Global *tsdata){
	/* update gamma_solver data - prepare new linear term */

	/* global data */
	Vec M_global = tsdata->get_thetavector()->get_vector(); /* parameters of the model */
	Vec x_global = tsdata->get_datavector()->get_vector(); /* original data */
	Vec b_global = gammadata->get_b()->get_vector(); /* rhs of QP gamma problem */

	/* local data */
	const double *theta;
	double *b_local_arr;

	/* get local arrays */
	TRY( VecGetArrayRead(M_global,&theta) );
	TRY( VecGetArray(b_global,&b_local_arr) );	

	/* get constants of the problem */
	int K = this->Klocal;
	int xdim = this->xdim;
	int xmem = this->xmemlocal;
	int xmem_max = max_array(GlobalManager.get_size(), this->xmem);
	int blocksize = 1 + xdim*xmem;
	int gammak_size = T-xmem;
	int theta_start; /* where in theta start actual coefficients */
	int k,t_mem,t,n,i;
	double Ax, value;
	double x_model[xdim];

	int t_scatter = 50; /* > xmem, how long time-series to scatter to all processors */
	int t_in_scatter;
	int is_begin = 0;
	int is_end = min(is_begin + t_scatter,T);

	VecScatter ctx;
	IS scatter_is;
	IS scatter_is_to;
	Vec x_scatter;
	TRY( VecCreateSeq(PETSC_COMM_SELF, t_scatter*xdim, &x_scatter) );
	const double *x_scatter_arr;
	while(is_end <= T && is_begin < is_end){ /* while there is something to scatter */
		/* scatter part of time serie */
		TRY( ISCreateStride(PETSC_COMM_WORLD, (is_end-is_begin)*xdim, is_begin*xdim, 1, &scatter_is) );
		TRY( ISCreateStride(PETSC_COMM_SELF, (is_end-is_begin)*xdim, 0, 1, &scatter_is_to) );

		TRY( VecScatterCreate(tsdata->get_datavector()->get_vector(),scatter_is, x_scatter,scatter_is_to,&ctx) );
		TRY( VecScatterBegin(ctx,tsdata->get_datavector()->get_vector(),x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
		TRY( VecScatterEnd(ctx,tsdata->get_datavector()->get_vector(),x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
		TRY( VecScatterDestroy(&ctx) );

		TRY( ISDestroy(&scatter_is_to) );
		TRY( ISDestroy(&scatter_is) );

		TRY( PetscBarrier(NULL) );
				
		/* write begin of time-serie */
		TRY( VecGetArrayRead(x_scatter, &x_scatter_arr) );
		for(t=is_begin+xmem_max;t<is_end;t++){
			t_in_scatter = t - is_begin;

			/* x_scatter_arr[t_in_scatter*xdim] is a start */
			for(k=0;k<K;k++){
				theta_start = k*blocksize*xdim; /* where in theta start actual coefficients */

				for(n = 0; n < xdim; n++){
					/* mu */
					x_model[n] = theta[theta_start + n*blocksize]; 
				
					/* A */
					for(t_mem = 1; t_mem <= xmem; t_mem++){
						/* add multiplication with A_{t_mem} */
						Ax = 0;
						for(i = 0; i < xdim; i++){
							Ax += theta[theta_start + n*blocksize + 1 + (t_mem-1)*xdim + i]*x_scatter_arr[(t_in_scatter-t_mem)*xdim+i]; 
						}
						x_model[n] += Ax;
					}
					
					/* B */
					//todo: process u(t)
				}

				/* now I have computed x_model(t), therefore I just compute || x_model(t)-x_data(t) ||^2 */
				value = 0;
				for(n=0;n<xdim;n++){
					value += (x_model[n] - x_scatter_arr[t_in_scatter*xdim+n])*(x_model[n] - x_scatter_arr[t_in_scatter*xdim+n]);
				}

				/* write the error into linear term */
				value = (-1)*value;
				b_local_arr[k*gammak_size+t-xmem] = value;
			}

		}
		TRY( VecRestoreArrayRead(x_scatter, &x_scatter_arr) );
				
		/* update upper and lower index of scatter_is */
		is_begin = is_begin + t_scatter - xmem_max;
		is_end = min(is_begin + t_scatter,T);

	}
	TRY( VecDestroy(&x_scatter) );

	TRY( VecRestoreArrayRead(M_global,&theta) );
	TRY( VecRestoreArray(b_global,&b_local_arr) );	

}


void VarxH1FEMModel_Global::scatter_xn(Vec &x_global, Vec &xn_local, int xdim, int n){
	/* compute T */
	int x_global_size;
	TRY( VecGetSize(x_global,&x_global_size) );
	int T = x_global_size/(double)xdim;

	VecScatter ctx; 
	IS xn_local_is;

	TRY( ISCreateStride(PETSC_COMM_WORLD, T, n, xdim, &xn_local_is) );

	TRY( VecScatterCreateToAll(x_global,&ctx,&xn_local) );
	TRY( VecScatterBegin(ctx, x_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
	TRY( VecScatterEnd(ctx, x_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
	TRY( VecScatterDestroy(&ctx) );

	TRY( ISDestroy(&xn_local_is) );
}

/* update theta solver */
void VarxH1FEMModel_Global::update_thetasolver(GeneralSolver *thetasolver, const TSData_Global *tsdata){
	/* update theta solver - prepare new BlockDiag matrix data and right-hand side vector */	

	/* pointers to global data */
	Vec x_global = tsdata->get_datavector()->get_vector();
	
	/* get constants of the problem */
	int Klocal = this->Klocal;
	int T = this->T;
	int xdim = this->xdim;
	int xmem = this->xmemlocal;
//	int umem = this->umem;
	int blocksize = 1 + xdim*xmem; /* see get_thetavectorlength() */

	/* blocks of the matrix */
	BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *A = dynamic_cast<BlockDiagMatrix<PetscVector,LocalDenseMatrix<PetscVector> > *>(thetadata->get_A());
	LocalDenseMatrix<PetscVector> **blocks = A->get_blocks();

	/* scatter and index sets with components of x corresponding to one dimension */
	Vec xn1_vec;  /* scattered xn_global to each processor */
	Vec xn2_vec;
	TRY( VecCreateSeq(PETSC_COMM_SELF, T, &xn1_vec) );
	TRY( VecCreateSeq(PETSC_COMM_SELF, T, &xn2_vec) );

	int gammak_size = T-xmem;

	/* prepare local gammas */
	Vec gamma_local;
	TRY( VecCreateSeq(PETSC_COMM_SELF, Klocal*gammak_size, &gamma_local) );
	TRY( VecGetLocalVector(gammadata->get_x()->get_vector(),gamma_local) );
	
	int Kmax = max_array(num, K);
	Vec gammak_vecs[Kmax];
	IS gammak_is[Kmax];
	int k;
	for(k=0;k<Kmax;k++){
		if(k < Klocal){
			TRY( ISCreateStride(PETSC_COMM_SELF, gammak_size, k*gammak_size, 1, &gammak_is[k]) );
		} else {
			TRY( ISCreateStride(PETSC_COMM_SELF, 0, k*gammak_size, 1, &gammak_is[k]) );
		}
		TRY( VecGetSubVector(gamma_local, gammak_is[k], &(gammak_vecs[k])) );
	}
	
	TRY( PetscBarrier(NULL) );
	
	/* prepare local b */
	double *b_arr;
	//TODO: temp
	TRY( VecSet(thetadata->get_b()->get_vector(), 0.0) );
	
	TRY( VecGetArray(thetadata->get_b()->get_vector(), &b_arr) );
		
	/* index set for manipulating with local parts of xn according to xmem */
	IS xn1sub_is;
	IS xn2sub_is;
	Vec xn1sub_vec; 
	Vec xn2sub_vec;
	int xn1sub_nmb;
	int xn2sub_nmb;
	
	/* vector with xn1sub.*gammak_vecs */
	Vec xn2subgammak_vec;
	TRY( VecCreateSeq(PETSC_COMM_SELF, gammak_size, &xn2subgammak_vec) );
		
	/* go through clusters and fill matrices */
	int col,row, xdim1,xdim2, i,j;
	double value;

	double coeff = 1;

	/* matrix */
	for(xdim1 = 0; xdim1<xdim+1; xdim1++){
//	for(xdim1 = 0; xdim1<1; xdim1++){
		/* constant */
		if(xdim1 == 0){
			/* the first row is 1 */
			TRY( VecSet(xn1_vec, coeff) );
			VecAssemblyBegin(xn1_vec);
			VecAssemblyEnd(xn1_vec);
			
			xn1sub_nmb = 1;
		}
		
		/* x rows */
		if(xdim1 >= 1 && xdim1 < 1+xdim){
			/* scatter xn with n=xdim1 */
			scatter_xn(x_global, xn1_vec, xdim, xdim1-1);
			xn1sub_nmb = xmem;

			TRY( VecView(x_global, PETSC_VIEWER_STDOUT_SELF) );
			TRY( VecView(xn1_vec, PETSC_VIEWER_STDOUT_SELF) );

		}

		TRY( PetscBarrier(NULL) );

		for(xdim2 = 0; xdim2<xdim+1;xdim2++){
//		for(xdim2 = 0; xdim2<1;xdim2++){

			/* constant */
			if(xdim2 == 0){
				/* the first row is 1 */
				TRY( VecSet(xn2_vec, coeff) );
				VecAssemblyBegin(xn2_vec);
				VecAssemblyEnd(xn2_vec);

				xn2sub_nmb = 1;
			}
		
			/* x rows */
			if(xdim2 >= 1 && xdim2 < 1+xdim){
				/* scatter xn with n=xdim2 */
				scatter_xn(x_global, xn2_vec, xdim, xdim2-1);
				xn2sub_nmb = xmem;
			}
			
			/* now I have two vectors with components corresponding with one dimension, what to do next? - local computations */
			
			coutMaster << "computing with: " << xdim1 << ", " << xdim2 << std::endl;
			
			for(i=0;i<xn1sub_nmb;i++){
				/* prepare subvector from xn_vec which corresponds to row of Z */
				TRY( ISCreateStride(PETSC_COMM_SELF, gammak_size, i, 1, &xn1sub_is) );
				TRY( VecGetSubVector(xn1_vec, xn1sub_is, &xn1sub_vec) );				
				
//				TRY( VecView(xn1_vec, PETSC_VIEWER_STDOUT_SELF) );
//				TRY( VecView(xn1sub_vec, PETSC_VIEWER_STDOUT_SELF) );
				
				/* ----- MATRIX --- */
				for(j=0;j<xn2sub_nmb;j++){
					/* prepare subvector from xn_vec which corresponds to column of Z */
					TRY( ISCreateStride(PETSC_COMM_SELF, gammak_size, j, 1, &xn2sub_is) );
					TRY( VecGetSubVector(xn2_vec, xn2sub_is, &xn2sub_vec) );				
					
					/* go through clusters and fill blocks */
					for(k=0;k<Klocal;k++){
						TRY( VecPointwiseMult(xn2subgammak_vec, xn2sub_vec, gammak_vecs[k]) );

						/* compute dot product xn1sub_vec*xn2sub_vec */
						TRY( VecDot(xn2subgammak_vec,xn1sub_vec,&value) );

						/* write values to symmetric matrix */
						if(xdim1 == 0){
							row = xdim1 + (xn1sub_nmb-i-1)*xdim;
						}
						if(xdim2 == 0){
							col = xdim2 + (xn2sub_nmb-j-1)*xdim;
						}
						
						if(xdim1 >= 1 && xdim1 < 1+xdim){
							row = xdim1 + (xn1sub_nmb-i-1)*xdim;
						}						
						if(xdim2 >= 1 && xdim2 < 1+xdim){
							col = xdim2 + (xn2sub_nmb-j-1)*xdim;
						}						


						if(row == col){
							/* diagonal entry */
							blocks[k*xdim]->set_value(row,col,value);
						} else {
							/* nondiagonal entry */
							blocks[k*xdim]->set_value(row,col,value);
//							blocks[k*xdim]->add_value(col,row,value);
						}
					}
					
					/* restore, clear and prepare for next row */
					TRY( VecRestoreSubVector(xn2_vec, xn2sub_is, &xn2sub_vec) );
					TRY( ISDestroy(&xn2sub_is) );

					coutMaster << " -i,j: " << i << ", " << j << std::endl;

				}
				
				/* ----- VECTOR ----- */
				if(xdim2 >= 1 && xdim2 < 1+xdim){
					TRY( ISCreateStride(PETSC_COMM_SELF, gammak_size, xmem, 1, &xn2sub_is) );
					TRY( VecGetSubVector(xn2_vec, xn2sub_is, &xn2sub_vec) );
					/* go through clusters and fill RHS vector */
					for(k=0;k<Klocal;k++){
						TRY( VecPointwiseMult(xn2subgammak_vec, xn2sub_vec, gammak_vecs[k]) );
//						MyVecPointwiseMult(xn1subgammak_vec, xn1sub_vec, gammak_vecs[k]);

						/* compute dot product xn1sub_vec*xn2sub_vec */
						TRY( VecDot(xn2subgammak_vec,xn1sub_vec,&value) );
//						MyVecDot(xn1subgammak_vec,xn2sub_vec,&value);

						if(xdim1 == 0){
							row = k*blocksize*xdim + (xdim2-1)*blocksize;
						}

						/* x rows */
						if(xdim1 >= 1 && xdim1 < 1+xdim){
							row = k*blocksize*xdim + (xdim2-1)*blocksize + xdim1 + (xn1sub_nmb-i-1)*xdim;

//							row = k*blocksize*xdim + (xn1sub_nmb-i-1)*blocksize + xdim2;
//							row = k*blocksize*xdim + i*xdim + (xdim2-1)*blocksize;
//							row = k*blocksize*xdim + i*blocksize + xdim2;

						}

						b_arr[row] = value;

					}
					TRY( VecRestoreSubVector(xn2_vec, xn2sub_is, &xn2sub_vec) );
					TRY( ISDestroy(&xn2sub_is) );
				}
				
				/* restore, clear and prepare for next row */
				TRY( VecRestoreSubVector(xn1_vec, xn1sub_is, &xn1sub_vec) );
				TRY( ISDestroy(&xn1sub_is) );
			}
			
			
			
		}
		
	}

	TRY( VecDestroy(&xn2subgammak_vec) );
	coutAll << "b = ";
	print_array(coutAll,Klocal*xdim*blocksize,b_arr);
	coutAll << std::endl;
	coutAll.synchronize();

	TRY( VecRestoreArray(thetadata->get_b()->get_vector(), &b_arr) );
//	TRY( VecRestoreLocalVector(gammadata->get_x()->get_vector(), gamma_local));

	/* restore gammas */
	for(k=0;k<Kmax;k++){
		TRY( VecRestoreSubVector(gammadata->get_x()->get_vector(), gammak_is[k], &(gammak_vecs[k])) );
		TRY( ISDestroy(&gammak_is[k]) );
	}

//	TRY( VecDestroy(&gamma_local) );
//	TRY( VecDestroy(&xn1_vec) );
//	TRY( VecDestroy(&xn2_vec) );

	/* assemble blocks of matrix */
	for(k=0;k<Klocal;k++){
		blocks[k*xdim]->assemble();
	}

	// TODO: temp print
	coutMaster << "------------- TEMP PRINT ------------" << std::endl;
	coutAll << "proc: " << GlobalManager.get_rank() << std::endl;
	for(i=0;i< Klocal*xdim;i++){
		coutMaster.push();
		coutAll << "block:" << i << std::endl;
		blocks[i]->printcontent(coutAll);
		coutMaster.pop();
	}
	coutAll.synchronize();
	
	/* regularization */
/*	for(k=0;k<Klocal;k++){
		for(i=0;i<blocksize;i++){
			
		}
	}
	*/
	/* rhs */
//	for(n=0;n<xdim;n++){
		/* get global xn_global */
//		TRY( ISCreateStride(PETSC_COMM_WORLD, T, n, xdim, &xn_is) );
//		TRY( VecGetSubVector(X,xn_is,&xn_global) );
		
		/* scatter xn_global to xn_local */
//		TRY( VecScatterCreateToAll(xn_global,&ctx,&xn_local) );
//		TRY( VecScatterBegin(ctx, xn_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
//		TRY( VecScatterEnd(ctx, xn_global, xn_local, INSERT_VALUES, SCATTER_FORWARD) );
//		TRY( VecScatterDestroy(&ctx) );

		/* now I have my own seq vector xn_local with x(n,:) */
//		TRY( VecGetArrayRead(xn_local, &xn_local_arr) );
		
//		for(k=0;k<Klocal;k++){
			/* RHS - compute dot product dot( X(xmem + n*(xmem+T), xmem+T-1  + n*(xmem+T)), gamma(k*T,(k+1)*T-1) ); */
//			value = 0;
//			for(t=0;t<T;t++){
//				value += xn_local_arr[t]*gamma_arr[k*T+t];
//			}
//			b_arr[k*xdim + n] = value;
//		}	

//		TRY( VecRestoreArrayRead(xn_local, &xn_local_arr) );

		/* restore subvector with xn_global */
//		TRY( VecRestoreSubVector(X,xn_is,&xn_global) );
//		TRY( ISDestroy(&xn_is) );

		
//	}



	
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
