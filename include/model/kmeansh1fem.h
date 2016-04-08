#ifndef PASC_KMEANSH1FEMMODEL_H
#define	PASC_KMEANSH1FEMMODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "model/tsmodel.h"
#include "matrix/blockdiaglaplace_explicit.h"

#include "feasibleset/simplex.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

#include "solver/diagsolver.h"
#include "data/diagdata.h"

#include "data/tsdata.h"

namespace pascinference {

/* KMEANSH1MODEL */ 
template<class VectorBase>
class KmeansH1FEMModel: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata;
		DiagData<VectorBase> *thetadata;

		double penalty;
	public:
		KmeansH1FEMModel(int T, int dim, int K, double penalty);
		~KmeansH1FEMModel();

		void print(std::ostream &output) const;
		std::string get_name() const;

		int get_datavectorlength();
		int get_gammavectorlength();
		int get_thetavectorlength();
		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);

		void update_gammasolver(GeneralSolver *gamma_solver, const TSData<VectorBase> *tsdata);
		void update_thetasolver(GeneralSolver *theta_solver, const TSData<VectorBase> *tsdata);
	
		void generate_data(TSData<VectorBase> *tsdata);

		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata);

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
KmeansH1FEMModel<VectorBase>::KmeansH1FEMModel(int newT, int newdim, int newK, double penalty) {
	if(DEBUG_MODE >= 100) coutMaster << "(KmeansH1FEMModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;

	this->penalty = penalty;
}

/* destructor */
template<class VectorBase>
KmeansH1FEMModel<VectorBase>::~KmeansH1FEMModel(){
	if(DEBUG_MODE >= 100) coutMaster << "(KmeansH1FEMModel)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about problem */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:    " << this->T << std::endl;
	output <<  " - K:    " << this->K << std::endl;
	output <<  " - dim:  " << this->dim << std::endl;
	output <<  " - penalty:  " << this->penalty << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string KmeansH1FEMModel<VectorBase>::get_name() const {
	return "KMEANS-H1-FEM Time-Series Model";	
}

/* get the appropriate length of datavector */
template<class VectorBase>
int KmeansH1FEMModel<VectorBase>::get_datavectorlength(){
	return this->dim*this->T;
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int KmeansH1FEMModel<VectorBase>::get_gammavectorlength(){
	return this->K*this->T;
}

/* get the appropriate length of thetavector */
template<class VectorBase>
int KmeansH1FEMModel<VectorBase>::get_thetavectorlength(){
	return this->K*this->dim;
}


/* prepare gamma solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<VectorBase>();
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<VectorBase>(*gammadata->get_x0())); /* create new linear term of QP problem */
	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<VectorBase>(*gammadata->get_x0(),this->K, this->penalty)); /* create new blockdiagonal matrix */
	gammadata->set_feasibleset(new SimplexFeasibleSet<VectorBase>(this->T,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<VectorBase>(*gammadata);
	
	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

}

/* prepare theta solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	/* in this case, theta problem is system with diagonal matrix */
	
	/* create data */
	thetadata = new DiagData<VectorBase>();
	thetadata->set_a(new GeneralVector<VectorBase>(*tsdata->get_thetavector()));
	thetadata->set_b(new GeneralVector<VectorBase>(*tsdata->get_thetavector()));
	thetadata->set_x(tsdata->get_thetavector());

	/* create solver */
	*thetasolver = new DiagSolver<VectorBase>(*thetadata);
	
}

/* destroy gamma solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
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
void KmeansH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	
	/* destroy data */
	free(thetadata->get_a());
	free(thetadata->get_b());
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

}

template<class VectorBase>
double KmeansH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

/* ---------------------- GENERAL ----------------------- */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::update_gammasolver(GeneralSolver *gammasolver, const TSData<VectorBase> *tsdata){
	/* update gamma_solver data - prepare new linear term */

	typedef GeneralVector<VectorBase> (&pVector);

	/* pointers to data */
	pVector theta = *(tsdata->get_thetavector());
	pVector data = *(tsdata->get_datavector());
	pVector b = *(gammadata->get_b());

	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	
	int t,n,k;
	double dot_sum, value;
	
	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			dot_sum = 0.0;
			for(n=0;n<dim;n++){
				value = data(n*T+t) - theta(k*dim + n); //TODO: problem with petscvector (double = subvector)
				dot_sum += value*value;
			}
			b(k*T + t) = -dot_sum;
		}
	}	
}

/* update theta solver */
template<class VectorBase>
void KmeansH1FEMModel<VectorBase>::update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	/* update theta solver - prepare new matrix vector and right-hand side vector */	
	
	typedef GeneralVector<VectorBase> (&pVector);

	/* pointers to data */
	pVector gamma = *(tsdata->get_gammavector());
	pVector data = *(tsdata->get_datavector());
	pVector a = *(thetadata->get_a());
	pVector b = *(thetadata->get_b());
	
	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	
	int k,i;
	double sum_gammak, gTd;
	for(k=0;k<K;k++){

		/* compute sum of gamma[k] */
		sum_gammak = sum(gamma(k*T,(k+1)*T-1));
		
		/* "a" has to be non-zero */
		if(sum_gammak == 0.0){
			sum_gammak = 1.0;
		}

		for(i=0;i<this->dim;i++){
			
			gTd = dot(gamma(k*T,(k+1)*T-1),data(i*T,(i+1)*T-1));

			a(k*dim+i) = sum_gammak;
			b(k*dim+i) = gTd;
		}
	}

}

/* ---------------------- PETSCVECTOR ----------------------- */
#ifdef USE_PETSCVECTOR

typedef petscvector::PetscVector GlobalPetscVector;

/* update gamma solver */ //TODO: inplement for each type
template<>
void KmeansH1FEMModel<GlobalPetscVector>::update_gammasolver(GeneralSolver *gammasolver, const TSData<GlobalPetscVector> *tsdata){
	/* update gamma_solver data - prepare new linear term */

	typedef GeneralVector<GlobalPetscVector> (&pVector);

	/* pointers to data */
	pVector theta = *(tsdata->get_thetavector());
	pVector data = *(tsdata->get_datavector());
	pVector b = *(gammadata->get_b());

	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	
	int t,k;
	
	Vec data_n;
	Vec theta_n;
	IS isdata_n, istheta_n;
	
	double temp_dot1,temp_dot2,temp_dot3;
	
	for(t=0;t<T;t++){
		TRY( ISCreateStride(PETSC_COMM_WORLD, dim, t, T, &isdata_n));
		TRY( VecGetSubVector(data.get_vector(),isdata_n, &data_n) );

		for(k=0;k<K;k++){
			TRY( ISCreateStride(PETSC_COMM_WORLD, dim, k*dim, 1, &istheta_n));
			TRY( VecGetSubVector(theta.get_vector(),istheta_n, &theta_n) );

			/* dot(data_n - theta_n,data_n - theta_n) */
			TRY( VecDot(data_n,data_n,&temp_dot1) );
			TRY( VecDot(data_n,theta_n,&temp_dot2) );
			TRY( VecDot(theta_n,theta_n,&temp_dot3) );

			b.set(k*T + t, -temp_dot1 + 2*temp_dot2 - temp_dot3);

			TRY( VecRestoreSubVector(theta.get_vector(),istheta_n, &theta_n) );
			TRY( ISDestroy(&istheta_n) );
		}
		
		TRY( VecRestoreSubVector(data.get_vector(),isdata_n, &data_n) );
		TRY( ISDestroy(&isdata_n) );
		
	}	
}

/* update theta solver */
template<>
void KmeansH1FEMModel<GlobalPetscVector>::update_thetasolver(GeneralSolver *thetasolver, const TSData<GlobalPetscVector> *tsdata){
	/* update theta solver - prepare new matrix vector and right-hand side vector */	
	
	typedef GeneralVector<GlobalPetscVector> (&pVector);

	/* pointers to data */
	pVector gamma = *(tsdata->get_gammavector());
	pVector data = *(tsdata->get_datavector());
	pVector a = *(thetadata->get_a());
	pVector b = *(thetadata->get_b());
	
	int K = this->K;
	int T = this->T;
	int dim = this->dim;
	
	int k,i;
	double sum_gammak, gTd;

	IS is, is2; 
	Vec temp1, temp2;

	for(k=0;k<K;k++){

		TRY( ISCreateStride(PETSC_COMM_WORLD, T, k*T, 1, &is));
		TRY( VecGetSubVector(gamma.get_vector(),is, &temp1) );

		/* compute sum of gamma[k] */
		TRY( VecSum(temp1, &sum_gammak) );
		
		/* "a" has to be non-zero */
		if(sum_gammak == 0.0){
			sum_gammak = 1.0;
		}

		for(i=0;i<this->dim;i++){

			TRY( ISCreateStride(PETSC_COMM_WORLD, T, i*T, 1, &is2));
			TRY( VecGetSubVector(data.get_vector(),is2, &temp2) );

			TRY( VecDot(temp1, temp2, &gTd) );

			a.set(k*dim+i, sum_gammak);
			b.set(k*dim+i, gTd);

			TRY( VecRestoreSubVector(data.get_vector(),is2, &temp2) );
			TRY( ISDestroy(&is2) );

		}

		TRY( VecRestoreSubVector(gamma.get_vector(),is, &temp1) );
		TRY( ISDestroy(&is) );

	}

}

#endif




} /* end namespace */

#endif
