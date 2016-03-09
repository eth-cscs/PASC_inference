#ifndef KMEANSH1MODEL_H
#define	KMEANSH1MODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "model/tsmodel.h"
#include "matrix/blockdiaglaplace_explicit.h"

#include "feasibleset/simplex.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"
#include "result/qpresult.h"

#include "data/tsdata.h"
#include "result/tsresult.h"

namespace pascinference {

/* KMEANSH1MODEL */ 
template<class VectorBase>
class KmeansH1Model: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata;
		QPResult<VectorBase> *gammaresult;

		QPData<VectorBase> *thetadata;
		QPResult<VectorBase> *thetaresult;

	public:
		KmeansH1Model(int T, int dim, int K);
		~KmeansH1Model();

		void print(std::ostream &output) const;
		std::string get_name() const;

		int get_datavectorlength();
		int get_gammavectorlength();
		int get_thetavectorlength();
		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult);
		
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
KmeansH1Model<VectorBase>::KmeansH1Model(int newT, int newdim, int newK) {
	if(DEBUG_MODE >= 100) std::cout << "(KmeansH1Model)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;

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
void KmeansH1Model<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult){
	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<VectorBase>();
	gammadata->x0 = tsresult->get_gammavector(); /* the initial approximation of QP problem is gammavector */
	gammadata->b = new GeneralVector<VectorBase>(*gammadata->x0); /* create new linear term of QP problem */
	gammadata->A = new BlockDiagLaplaceExplicitMatrix<VectorBase>(*gammadata->x0,this->K); /* create new blockdiagonal matrix */
	gammadata->feasibleset = new SimplexFeasibleSet<VectorBase>(this->T,this->K); /* the feasible set of QP is simplex */ 	

	/* create results */
	gammaresult = new QPResult<VectorBase>();
	gammaresult->x = tsresult->get_gammavector(); /* the solution of QP problem is gamma */

	/* create solver */
	*gammasolver = new QPSolver<VectorBase>(*gammadata,*gammaresult);
	
}

/* prepare theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult){
	/* in this case, theta problem is QP with empty feasible set */
	// TODO: write this funny function
	
	/* create data */
	thetadata = new QPData<VectorBase>();

	/* create results */
	thetaresult = new QPResult<VectorBase>();

	/* create solver */
	*thetasolver = NULL;
//	*thetasolver = new QPSolver<VectorBase>(*gammadata,*gammaresult);
	
}

/* destroy gamma solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult){
	/* I created this objects, I should destroy them */

	/* destroy data */
	free(gammadata->b);
	free(gammadata->A);
	free(gammadata->feasibleset);
	free(gammadata);

	/* destroy results */
	free(gammaresult);

	/* destroy solver */
	free(*gammasolver);

	
}

/* destroy theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata, const TSResult<VectorBase> *tsresult){
	/* in this case, theta problem is QP with empty feasible set */
	
	/* destroy data */
	free(thetadata);

	/* destroy results */
	free(thetaresult);

	/* destroy solver */
//	free(*thetasolver);

	
}



} /* end namespace */

#endif
