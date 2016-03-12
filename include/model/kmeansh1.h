#ifndef PASC_KMEANSH1MODEL_H
#define	PASC_KMEANSH1MODEL_H

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
class KmeansH1Model: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata;
		DiagData<VectorBase> *thetadata;

		double penalty;
	public:
		KmeansH1Model(int T, int dim, int K, double penalty);
		~KmeansH1Model();

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
KmeansH1Model<VectorBase>::KmeansH1Model(int newT, int newdim, int newK, double penalty) {
	if(DEBUG_MODE >= 100) std::cout << "(KmeansH1Model)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;

	this->penalty = penalty;
}

/* destructor */
template<class VectorBase>
KmeansH1Model<VectorBase>::~KmeansH1Model(){
	if(DEBUG_MODE >= 100) std::cout << "(KmeansH1Model)DESTRUCTOR" << std::endl;
	
	/* destroy auxiliary vectors */

}


/* print info about problem */
template<class VectorBase>
void KmeansH1Model<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:    " << this->T << std::endl;
	output << "  - K:    " << this->K << std::endl;
	output << "  - dim:  " << this->dim << std::endl;
	output << "  - penalty:  " << this->penalty << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string KmeansH1Model<VectorBase>::get_name() const {
	return "KMEANS-H1 Time-Series Model";	
}

/* get the appropriate length of datavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_datavectorlength(){
	return this->dim*this->T;
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_gammavectorlength(){
	return this->K*this->T;
}

/* get the appropriate length of thetavector */
template<class VectorBase>
int KmeansH1Model<VectorBase>::get_thetavectorlength(){
	return this->K*this->dim;
}


/* prepare gamma solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
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
//	gammadata->get_x0()->set_random();
	int clusterT = ceil((double)this->T/(double)this->K);
	
	/* generate random data */
	int k,t;
	for(k=0;k<this->K;k++){ /* go through clusters */ 
		for(t=k*clusterT;t < std::min((k+1)*clusterT,this->T);t++){ /* go through part of time serie corresponding to this cluster */
		
			if(k==0){
				(*(tsdata->get_gammavector()))(t) = 1.0;
			}
			if(k==1){
				(*(tsdata->get_gammavector()))(t+this->T) = 1.0;
			}
			if(k==2){
				(*(tsdata->get_gammavector()))(t+2*this->T) = 1.0;
			}
			
		}
	}
	
	
	/* project random values to feasible set */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	std::cout << *(gammadata->get_x0()) << std::endl;



}

/* prepare theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
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
void KmeansH1Model<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
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
void KmeansH1Model<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	
	/* destroy data */
	free(thetadata->get_a());
	free(thetadata->get_b());
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

}

/* update gamma solver */ //TODO: inplement for each type
template<class VectorBase>
void KmeansH1Model<VectorBase>::update_gammasolver(GeneralSolver *gammasolver, const TSData<VectorBase> *tsdata){
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
		for(t=0;t<T;t++){ // TODO: this could be performed parallely 
			dot_sum = 0.0;
			for(n=0;n<dim;n++){
				value = data.get(n*T+t) - theta.get(k*dim + n);
				dot_sum += value*value;
			}
			b.set(k*T + t, dot_sum);
//			b(k*T + t) = dot_sum;
		}
	}	
}

/* update theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
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

//		std::cout << "gamma_" << k << ":" << gamma(k*T,(k+1)*T-1) << std::endl;

		for(i=0;i<this->dim;i++){
//			std::cout << "data_" << i << ":" << data(i*T,(i+1)*T-1) << std::endl;

			
			gTd = dot(gamma(k*T,(k+1)*T-1),data(i*T,(i+1)*T-1));

			a.set(k*dim+i, sum_gammak);
			b.set(k*dim+i, gTd);
		}
	}

}

template<class VectorBase>
double KmeansH1Model<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}


} /* end namespace */

#endif
