#ifndef PASC_KMEANSH1FEMMODEL_H
#define	PASC_KMEANSH1FEMMODEL_H

#ifndef USE_PETSCVECTOR
 #error 'KMEANSH1FEMMODEL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

/* gamma problem */
#include "matrix/blocklaplace.h"
//#include "matrix/blockdiaglaplace_vector.h"

#include "feasibleset/simplex_local.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

/* theta problem */
#include "solver/diagsolver.h"
#include "data/diagdata.h"

#include "data/kmeansdata.h"


namespace pascinference {

template<class VectorBase>
class KmeansH1FEMModel: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex, will be solved by SPG-QP */
	 	DiagData<VectorBase> *thetadata; /**< QP with diagonal matrix of dimension xdim*K, will be solved by diagsolver  */

		/* model specific variables */
		double epssqr; /**< penalty coeficient */
		
	public:
	
		KmeansH1FEMModel(KmeansData<VectorBase> &tsdata, int xdim, int K, double epssqr);
		~KmeansH1FEMModel();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		std::string get_name() const;
		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);

		void update_gammasolver(GeneralSolver *gamma_solver, const TSData<VectorBase> *tsdata);
		void update_thetasolver(GeneralSolver *theta_solver, const TSData<VectorBase> *tsdata);
	
		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata);

		QPData<VectorBase> *get_gammadata() const;
		DiagData<VectorBase> *get_thetadata() const;

};

} // end of namespace


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<>
KmeansH1FEMModel<PetscVector>::KmeansH1FEMModel(KmeansData<PetscVector> &tsdata, int xdim, int K, double epssqr) {
	LOG_FUNC_BEGIN

	/* get original PETSc Vec from data vector */
	Vec data_Vec = tsdata.get_datavector()->get_vector();

	/* get T from data */
	TRY( VecGetSize(data_Vec, &(this->datavectorlength_global)) );
	this->T = this->datavectorlength_global/(double)xdim;

	/* get Tlocal from data */
	TRY( VecGetLocalSize(data_Vec, &(this->datavectorlength_local)) );
	this->Tlocal = this->datavectorlength_local/(double)xdim;

	/* get T ranges from data */
	int low, high;
	TRY( VecGetOwnershipRange(data_Vec, &low, &high) );
	this->Tbegin = low/(double)xdim;
	this->Tend = high/(double)xdim;

	/* set given parameters */
	this->xdim = xdim;
	this->K = K;
	this->epssqr = epssqr;

	int my_rank = GlobalManager.get_rank();

	/* compute vector lengths */

	/* prepage global vectors with gamma */
	this->gammavectorlength_local = K*this->Tlocal;
	this->gammavectorlength_global = K*this->T;

	/* prepare sequential vector with Theta - yes, all procesors will have the same information */
	this->thetavectorlength_local = xdim*K; 
	this->thetavectorlength_global = GlobalManager.get_size()*(this->thetavectorlength_local);

	/* set this model to data */
	tsdata.set_model(*this);
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
KmeansH1FEMModel<VectorBase>::~KmeansH1FEMModel(){
	LOG_FUNC_BEGIN
	
	/* destroy auxiliary vectors */
	
	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K:       " << this->K << std::endl;
	output <<  " - epssqr:  " << this->epssqr << std::endl;

	output <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;
	
	output.synchronize();	

	LOG_FUNC_END
}

/* print info about model */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info:  " << std::endl;
	output_global <<  "  - T:       " << this->T << std::endl;
	output_global <<  "  - xdim:    " << this->xdim << std::endl;

	output_global <<  " - K:       " << this->K << std::endl;
	output_global <<  " - epssqr:  " << this->epssqr << std::endl;

	output_global <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output_global <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output_global <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables:  " << std::endl;
	output_global.push();
	output_local << "Tlocal =" << std::setw(6) << this->Tlocal << " (" << this->Tbegin << "," << this->Tend << "), ";
	output_local << "datalength=" << std::setw(6) << this->datavectorlength_local << ", ";
	output_local << "gammalength=" << std::setw(6) << this->gammavectorlength_local << ", ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

/* print model solution */
template<>
void KmeansH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  "theta:" << std::endl;

	std::ostringstream temp;
	
	double *theta;
	TRY( VecGetArray(thetadata->get_x()->get_vector(),&theta) );
	
	int k,n;

	output_local.push();
	output_local << "- proc: " << GlobalManager.get_rank() << std::endl;
	output_local.push();
	for(k=0;k<this->K;k++){
		output_local <<  "- k = " << k << std::endl;

		/* mu */
		output_local.push();
		output_local <<  "- mu = [";
		for(n=0;n<xdim;n++){
			temp << theta[k*xdim + n];
			output_local << temp.str();
			if(n < xdim-1){
				output_local << ", ";
			}
			temp.str("");
		}
		output_local <<  "]" << std::endl;
		output_local.pop();

	}
	output_local.pop();
	output_local.pop();

	output_local.synchronize();
	output_global.synchronize();

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string KmeansH1FEMModel<VectorBase>::get_name() const {
	return "Kmeans-H1-FEM Time-Series Model";	
}

/* prepare gamma solver */
template<>
void KmeansH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN

	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<PetscVector>();
	
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */

	gammadata->set_A(new BlockLaplaceMatrix<PetscVector>(*(gammadata->get_x0()), this->K, this->epssqr)); 
	gammadata->set_feasibleset(new SimplexFeasibleSet_Local(this->Tlocal,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<PetscVector>(*gammadata);

	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	LOG_FUNC_END
}

/* prepare theta solver */
template<>
void KmeansH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN
	
	/* create data */
	thetadata = new DiagData<PetscVector>();
	thetadata->set_x(tsdata->get_thetavector()); /* the solution of QP problem is gamma */
	thetadata->set_a(new GeneralVector<PetscVector>(*thetadata->get_x())); /* create new diagonal of matrix */
	thetadata->set_b(new GeneralVector<PetscVector>(*thetadata->get_x())); /* create new rhs */

	/* create solver */
	*thetasolver = new DiagSolver<PetscVector>(*thetadata);
	
	LOG_FUNC_END
}

/* destroy gamma solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(gammadata->get_b());
	free(gammadata->get_A());
	free(gammadata->get_feasibleset());
	free(gammadata);

	/* destroy solver */
	free(*gammasolver);
	
	LOG_FUNC_END
}

/* destroy theta solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(thetadata->get_a());
	free(thetadata->get_b());
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

	LOG_FUNC_END
}

template<class VectorBase>
double KmeansH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

template<class VectorBase>
QPData<VectorBase>* KmeansH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
DiagData<VectorBase>* KmeansH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}

template<>
void KmeansH1FEMModel<PetscVector>::update_gammasolver(GeneralSolver *gammasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN

	/* update gamma_solver data - prepare new linear term */
	const double *theta_arr;
	TRY( VecGetArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	
	const double *data_arr;
	TRY( VecGetArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	
	double *b_arr;
	TRY( VecGetArray(gammadata->get_b()->get_vector(), &b_arr) );

	int k,t,n;
	double mydot;
	for(t=0;t<Tlocal;t++){
		for(k=0;k<K;k++){
			mydot = 0;
			for(n=0;n<xdim;n++){
				mydot += (data_arr[t*xdim+n] - theta_arr[k*xdim+n])*(data_arr[t*xdim+n] - theta_arr[k*xdim+n]); 
			}
			b_arr[k*Tlocal+t] = -mydot;
		}
	}

	/* restore arrays */
	TRY( VecRestoreArray(gammadata->get_b()->get_vector(), &b_arr) );
	TRY( VecRestoreArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	TRY( VecRestoreArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	
	LOG_FUNC_END
}

/* update theta solver */
template<>
void KmeansH1FEMModel<PetscVector>::update_thetasolver(GeneralSolver *thetasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN

	Vec gamma_Vec = tsdata->get_gammavector()->get_vector();
	Vec data_Vec = tsdata->get_datavector()->get_vector();

	/* get arrays from problem objects */
	double *a_arr;
	double *b_arr;
	TRY( VecGetArray(thetadata->get_a()->get_vector(), &a_arr) );
	TRY( VecGetArray(thetadata->get_b()->get_vector(), &b_arr) );

	Vec gammak_Vec;
	IS gammak_is;

	Vec datan_Vec;
	IS datan_is;
	
	int k,n;	
	double sum_gammak, gammaTx;
	/* through clusters */
	for(k=0;k<K;k++){ 
		/* get gammak */
		TRY( ISCreateStride(PETSC_COMM_WORLD, Tlocal, Tbegin*K + k*Tlocal, 1, &gammak_is) );
		TRY( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		/* compute sum of gammak */
		TRY( VecSum( gammak_Vec, &sum_gammak) );
		myround(sum_gammak,&sum_gammak);

		/* through dimensions */
		for(n=0;n<xdim;n++){
			/* get datan */
			TRY( ISCreateStride(PETSC_COMM_WORLD, Tlocal, Tbegin*xdim + n, xdim, &datan_is) );
			TRY( VecGetSubVector(data_Vec, datan_is, &datan_Vec) );
			
			/* compute dot product */
			TRY( VecDot(datan_Vec,gammak_Vec, &gammaTx) );
			myround(gammaTx, &gammaTx);
			
			/* store values to problem objects */
			a_arr[k*xdim+n] = sum_gammak;
			b_arr[k*xdim+n] = gammaTx;

			/* restore datan */
			TRY( VecRestoreSubVector(data_Vec, datan_is, &datan_Vec) );
			TRY( ISDestroy(&datan_is) );
		}

		/* restore gammak */
		TRY( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		TRY( ISDestroy(&gammak_is) );

	}	

	/* restore arrays to problem objects */
	TRY( VecRestoreArray(thetadata->get_a()->get_vector(), &a_arr) );
	TRY( VecRestoreArray(thetadata->get_b()->get_vector(), &b_arr) );

	LOG_FUNC_END
}



} /* end namespace */

#endif
