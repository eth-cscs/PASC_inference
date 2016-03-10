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

	public:
		KmeansH1Model(int T, int dim, int K);
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
void KmeansH1Model<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<VectorBase>();
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<VectorBase>(*gammadata->get_x0())); /* create new linear term of QP problem */
	gammadata->set_A(new BlockDiagLaplaceExplicitMatrix<VectorBase>(*gammadata->get_x0(),this->K)); /* create new blockdiagonal matrix */
	gammadata->set_feasibleset(new SimplexFeasibleSet<VectorBase>(this->T,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<VectorBase>(*gammadata);
	
}

/* prepare theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	/* in this case, theta problem is system with diagonal matrix */
	
	/* create data */
	thetadata = new DiagData<VectorBase>();
	thetadata->set_a(new GeneralVector<VectorBase>(*tsdata->get_thetavector()));
	thetadata->set_b(new GeneralVector<VectorBase>(*tsdata->get_thetavector()));

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

/* update gamma solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::update_gammasolver(GeneralSolver *gammasolver, const TSData<VectorBase> *tsdata){
	/* update gamma_solver data - prepare new linear term */

	
}

/* update theta solver */
template<class VectorBase>
void KmeansH1Model<VectorBase>::update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	/* update theta solver - prepare new matrix vector and right-hand side vector */	

}



} /* end namespace */

#endif
