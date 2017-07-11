#ifndef PASC_GRAPHH1FEMMODEL_H
#define PASC_GRAPHH1FEMMODEL_H

#include "general/common/common.h"
#include "general/algebra/fem/fem.h"

/* gamma problem */
//#include "algebra/matrix/blockgraphfree.h" // TODO: implement?
#include "general/algebra/matrix/blockgraphsparse.h"
#include "general/algebra/feasibleset/simplex_local.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"
#include "general/solver/spgqpsolver.h"
#include "general/solver/spgqpsolver_coeff.h"
#include "general/solver/permonsolver.h"
//#include "general/solver/taosolver.h"
#include "general/data/qpdata.h"

/* theta problem */
#include "general/solver/simplesolver.h"
#include "general/data/simpledata.h"
#include "general/data/tsdata.h"


/* default type of matrix: 0=FREE, 1=SPARSE */
#define GRAPHH1FEMMODEL_DEFAULT_MATRIX_TYPE 1
#define GRAPHH1FEMMODEL_DEFAULT_SCALEF 1

namespace pascinference {
namespace model {

/** \class GraphH1FEMModel
 *  \brief time-series model with quadratic penalty time-space regularisation.
 *
*/
template<class VectorBase>
class GraphH1FEMModel: public TSModel<VectorBase> {
	public:
		/** @brief type of solver used to solve inner gamma problem 
		 * 
		 * The inner problem leads to QP with equality constraints and bound constraints
		 * 
		 */
		typedef enum { 
			SOLVER_AUTO=0, /**< choose automatic solver */
			SOLVER_SPGQP=1, /**< CPU/GPU implementation of Spectral Projected Gradient method */
			SOLVER_SPGQP_COEFF=2, /**< CPU/GPU implementation of Spectral Projected Gradient method with special coefficient threatment */
			SOLVER_PERMON=3, /**< PERMONQP solver (augumented lagrangians combined with active-set method) */
			SOLVER_TAO=4 /**< TAO solver */
		} GammaSolverType;

	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex, will be solved by SPG-QP */
	 	SimpleData<VectorBase> *thetadata; /**< this problem is solved during assembly  */

		double epssqr; /**< penalty coeficient */
		
		/* for theta problem */
		GeneralMatrix<VectorBase> *A_shared; /**< matrix shared by gamma and theta solver */
		GeneralVector<VectorBase> *Agamma; /**< temp vector for storing A_shared*gamma */
		GeneralVector<VectorBase> *residuum; /**< temp vector for residuum computation */
		
		bool usethetainpenalty; /**< use the value of Theta in penalty parameter to scale blocks */
		bool scalef;			/**< divide whole function by T */
		
		GammaSolverType gammasolvertype; /**< the type of used solver */
		
		Fem<VectorBase> *fem;

	public:

		/** @brief constructor from data and regularisation constant
		 * 
		 * @param tsdata time-series data on which model operates
		 * @param epssqr regularisation constant
		 */ 	
		GraphH1FEMModel(TSData<VectorBase> &tsdata, double epssqr, Fem<VectorBase> *new_fem = NULL, bool usethetainpenalty = false);

		/** @brief destructor 
		 */ 
		~GraphH1FEMModel();

		/** @brief print info about model
		 * 
		 * @param output where to print
		 */	
		void print(ConsoleOutput &output) const;

		/** @brief print info about model
		 * 
		 * @param output_global where to print global info
		 * @param output_local where to print local info
		 */	
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief print solution of the model
		 * 
		 * @param output_global where to print global part
		 * @param output_local where to print local part
		 */	
		void printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		std::string get_name() const;
		
		void initialize_gammasolver(GeneralSolver **gamma_solver);
		void updatebeforesolve_gammasolver(GeneralSolver *gamma_solver);
		void updateaftersolve_gammasolver(GeneralSolver *gamma_solver);
		void finalize_gammasolver(GeneralSolver **gamma_solver);

		void initialize_thetasolver(GeneralSolver **theta_solver);
		void updatebeforesolve_thetasolver(GeneralSolver *theta_solver);
		void updateaftersolve_thetasolver(GeneralSolver *theta_solver);
		void finalize_thetasolver(GeneralSolver **theta_solver);
	
		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver);
		void get_linear_quadratic(double *linearL, double *quadraticL, GeneralSolver *gammasolver, GeneralSolver *thetasolver);
		
		QPData<VectorBase> *get_gammadata() const;
		SimpleData<VectorBase> *get_thetadata() const;
		BGMGraph<VectorBase> *get_graph() const;

		GeneralVector<VectorBase> *get_coordinatesVTK() const;
		int get_coordinatesVTK_dim() const;

		double get_aic(double L) const;
		
		void set_epssqr(double epssqr);
		bool get_usethetainpenalty() const;
		
		int get_T_reduced() const;
		int get_T() const;
		Decomposition<VectorBase> *get_decomposition_reduced() const;
		double get_fem_reduce() const;

};


}
} /* end of namespace */


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace model {

/* constructor */
template<class VectorBase>
GraphH1FEMModel<VectorBase>::GraphH1FEMModel(TSData<VectorBase> &new_tsdata, double epssqr, Fem<VectorBase> *new_fem, bool usethetainpenalty) {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
GraphH1FEMModel<VectorBase>::~GraphH1FEMModel(){
	LOG_FUNC_BEGIN
	
	/* destroy auxiliary vectors */
//	delete this->decomposition_reduced;
	
	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T                 : " << this->tsdata->get_T() << std::endl;
	output <<  " - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output <<  " - scalef            : " << printbool(this->scalef) << std::endl;

	/* information of reduced problem */
	output <<  " - fem_reduce        : " << this->fem->get_fem_reduce() << std::endl;
	output.push();
	this->fem->print(output,output);
	output.pop();
	
	output <<  " - K                 : " << this->tsdata->get_K() << std::endl;
	output <<  " - R                 : " << this->tsdata->get_R() << std::endl;
	output <<  " - epssqr            : " << this->epssqr << std::endl;
	output <<  " - usethetainpenalty : " << printbool(this->usethetainpenalty) << std::endl; 

	output <<  " - Graph             : " << std::endl;
	output.push();
	this->get_graph()->print(output);
	output.pop();

	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	output <<  " - gammasolvertype   : ";
	switch(this->gammasolvertype){
		case(SOLVER_AUTO): output << "AUTO"; break;
		case(SOLVER_SPGQP): output << "SPG-QP solver"; break;
		case(SOLVER_SPGQP_COEFF): output << "SPG-QP-COEFF solver"; break;
		case(SOLVER_PERMON): output << "PERMON QP solver"; break;
		case(SOLVER_TAO): output << "TAO QP solver"; break;
	}
	output << std::endl;
	
	output.synchronize();	

	LOG_FUNC_END
}

/* print info about model */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info" << std::endl;
	output_global <<  "  - T                 : " << this->tsdata->get_T() << std::endl;
	output_global <<  "  - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output_global <<  "  - scalef            : " << printbool(this->scalef) << std::endl;

	/* information of reduced problem */
	output_global <<  "  - fem_reduce        : " << this->fem->get_fem_reduce() << std::endl;
	output_global.push();
	this->fem->print(output_global, output_local);
	output_global.pop();

	output_global <<  "  - K                 : " << this->tsdata->get_K() << std::endl;
	output_global <<  "  - R                 : " << this->tsdata->get_R() << std::endl;
	output_global <<  "  - epssqr            : " << this->epssqr << std::endl;
	output_global <<  "  - usethetainpenalty : " << printbool(this->usethetainpenalty) << std::endl; 
	output_global <<  "  - gammasolvertype   : ";
	switch(this->gammasolvertype){
		case(SOLVER_AUTO): output_global << "AUTO"; break;
		case(SOLVER_SPGQP): output_global << "SPG-QP solver"; break;
		case(SOLVER_SPGQP_COEFF): output_global << "SPG-QP-COEFF solver"; break;
		case(SOLVER_PERMON): output_global << "PERMON QP solver"; break;
		case(SOLVER_TAO): output_global << "TAO QP solver"; break;
	}
	output_global << std::endl;

	output_global.push();
	this->get_graph()->print(output_global);
	output_global.pop();

	output_global <<  "  - thetalength       : " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables" << std::endl;
	output_global.push();
	output_local << "Tlocal =" << std::setw(6) << this->tsdata->get_Tlocal() << " (" << this->tsdata->get_Tbegin() << "," << this->tsdata->get_Tend() << "), ";
	output_local << "Rlocal =" << std::setw(6) << this->tsdata->get_Rlocal() << " (" << this->tsdata->get_Rbegin() << "," << this->tsdata->get_Rend() << "), ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

/* print model solution */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	//TODO

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string GraphH1FEMModel<VectorBase>::get_name() const {
	return "Graph-H1-FEM Time-Series Model";
}

/* set new penalty */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::set_epssqr(double epssqr) {
	LOG_FUNC_BEGIN

	this->epssqr = epssqr;

	/* use old T to scale the function to obtain the same scale of function values (idea from Olga) */
	double coeff = this->epssqr;
	if(this->scalef){
		coeff *= (1.0/((double)(this->get_T())));
	} else {
		coeff *= ((double)(this->get_T_reduced())/((double)(this->get_T())));
	}

	if(this->A_shared != NULL){
		/* SPARSE */
		((BlockGraphSparseMatrix<VectorBase>*)A_shared)->set_coeff(coeff);
	}

	LOG_FUNC_END
}

/* prepare gamma solver */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* prepare theta solver */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

/* destroy gamma solver */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */
	if(fem->is_reduced()){
		free(gammadata->get_x());
	}

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
void GraphH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(Agamma);
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

	LOG_FUNC_END
}

template<class VectorBase>
double GraphH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver){
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

template<class VectorBase> //TODO: ?
void GraphH1FEMModel<VectorBase>::get_linear_quadratic(double *linearL, double *quadraticL, GeneralSolver *gammasolver, GeneralSolver *thetasolver){
	*linearL = 0;
	*quadraticL = 0;
}

template<class VectorBase>
QPData<VectorBase>* GraphH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
SimpleData<VectorBase>* GraphH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}

template<class VectorBase>
BGMGraph<VectorBase> *GraphH1FEMModel<VectorBase>::get_graph() const {
	return this->tsdata->get_decomposition()->get_graph();
}

template<class VectorBase>
void GraphH1FEMModel<VectorBase>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void GraphH1FEMModel<VectorBase>::updateaftersolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


/* update theta solver */
template<class VectorBase>
void GraphH1FEMModel<VectorBase>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void GraphH1FEMModel<VectorBase>::updateaftersolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


template<class VectorBase>
GeneralVector<VectorBase> *GraphH1FEMModel<VectorBase>::get_coordinatesVTK() const{
	return this->get_graph()->get_coordinates();
}

template<class VectorBase>
int GraphH1FEMModel<VectorBase>::get_coordinatesVTK_dim() const{
	return this->get_graph()->get_dim();
}

template<class VectorBase>
double GraphH1FEMModel<VectorBase>::get_aic(double L) const{
	return 2*log(L) + this->tsdata->get_K();

}

template<class VectorBase>
bool GraphH1FEMModel<VectorBase>::get_usethetainpenalty() const{
	return this->usethetainpenalty;

}

template<class VectorBase>
int GraphH1FEMModel<VectorBase>::get_T_reduced() const {
	return fem->get_decomposition_reduced()->get_T();
}

template<class VectorBase>
int GraphH1FEMModel<VectorBase>::get_T() const {
	return this->tsdata->get_T();
}

template<class VectorBase>
double GraphH1FEMModel<VectorBase>::get_fem_reduce() const {
	return this->fem_reduce;
}

template<class VectorBase>
Decomposition<VectorBase> *GraphH1FEMModel<VectorBase>::get_decomposition_reduced() const {
	return this->fem->get_decomposition_reduced();
}



}
} /* end namespace */

#endif
