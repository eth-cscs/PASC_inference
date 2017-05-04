#ifndef PASC_ENTROPYH1FEMMODEL_H
#define PASC_ENTROPYH1FEMMODEL_H

#include "general/model/tsmodel.h"

/* gamma problem */
#include "general/algebra/matrix/blockgraphsparse.h"
#include "general/data/qpdata.h"
#include "general/algebra/feasibleset/simplex_local.h"
#include "general/solver/spgqpsolver.h"
#include "general/solver/spgqpsolver_coeff.h"
//#include "general/solver/taosolver.h"

#ifdef USE_PERMON
	#include "general/algebra/feasibleset/simplex_lineqbound.h"
	#include "general/solver/permonsolver.h"
#endif

/* theta problem */
#include "general/data/entropydata.h"
#include "general/solver/entropysolverdlib.h"
#include "general/solver/entropysolverspg.h"
#include "general/solver/entropysolvernewton.h"
#include "general/data/tsdata.h"

/* default type of matrix: 0=FREE, 1=SPARSE */
#define GRAPHH1FEMMODEL_DEFAULT_MATRIX_TYPE 1
#define GRAPHH1FEMMODEL_DEFAULT_SCALEF 1

namespace pascinference {
namespace model {

/** \class EntropyH1FEMModel
 *  \brief time-series model with quadratic penalty time-space regularisation and Entropy-based Theta model.
 *
*/
template<class VectorBase>
class EntropyH1FEMModel: public TSModel<VectorBase> {
	public:
		/** @brief type of solver used to solve inner gamma problem 
		 * 
		 * The inner problem leads to QP with equality constraints and bound constraints
		 * 
		 */
		typedef enum { 
			GSOLVER_AUTO=0,				/**< choose automatic solver */
			GSOLVER_SPGQP=1,			/**< CPU/GPU implementation of Spectral Projected Gradient method */
			GSOLVER_SPGQP_COEFF=2,		/**< CPU/GPU implementation of Spectral Projected Gradient method with special coefficient threatment */
			GSOLVER_PERMON=3,			/**< PERMONQP solver (augumented lagrangians combined with active-set method) */
			GSOLVER_TAO=4				/**< TAO solver */
		} GammaSolverType;

		/** @brief return name of gamma solver in string format
		 */
		std::string print_gammasolvertype(GammaSolverType gammasolvertype_in) const;

		/** @brief type of solver used to solve inner theta problem 
		 * 
		 * Maximum entropy - nonlinear optimization problem / system of nonlinear equations with integrals
		 * 
		 */
		typedef enum { 
			TSOLVER_AUTO=0,				/**< choose automatic solver */
			TSOLVER_ENTROPY_DLIB=1,		/**< use Dlib to solve integral problem */
			TSOLVER_ENTROPY_SPG=2,		/**< SPG (unconstrained) for solving integral problem */
			TSOLVER_ENTROPY_NEWTON=3	/**< Newton method for solving integral problem */
		} ThetaSolverType;

		/** @brief return name of theta solver in string format
		 */
		std::string print_thetasolvertype(ThetaSolverType thetasolvertype_in) const;

	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex, will be solved by SPG-QP */
	 	EntropyData<VectorBase> *thetadata; /**< for computing lambda-problem with integrals (Anna knows...)  */

		double epssqr; /**< penalty coeficient */
		int Km;			/**< number of moments */
		
		/* for theta problem */
		GeneralMatrix<VectorBase> *A_shared; /**< matrix shared by gamma and theta solver */
		GeneralVector<VectorBase> *residuum; /**< temp vector for residuum computation */
		
		bool scalef;			/**< divide whole function by T */
		
		GammaSolverType gammasolvertype; /**< the type of used solver for gamma problem */
		ThetaSolverType thetasolvertype; /**< the type of used solver for theta problem */
		
	public:

		/** @brief constructor from data and regularisation constant
		 * 
		 * @param tsdata time-series data on which model operates
		 * @param epssqr regularisation constant
		 */ 	
		EntropyH1FEMModel(TSData<VectorBase> &tsdata, int Km, double epssqr);

		/** @brief destructor 
		 */ 
		~EntropyH1FEMModel();

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
		EntropyData<VectorBase> *get_thetadata() const;
		BGMGraph<VectorBase> *get_graph() const;

		double get_aic(double L) const;
		int get_Km() const;
		void set_epssqr(double epssqr);
		int get_T() const;

};


}
} /* end of namespace */


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace model {

template<class VectorBase>
std::string EntropyH1FEMModel<VectorBase>::print_gammasolvertype(GammaSolverType gammasolvertype_in) const {
	std::string return_value = "?";
	switch(gammasolvertype_in){
		case(GSOLVER_AUTO): 			return_value = "AUTO"; break;
		case(GSOLVER_SPGQP): 			return_value = "SPG-QP solver"; break;
		case(GSOLVER_SPGQP_COEFF):	 	return_value = "SPG-QP-COEFF solver"; break;
		case(GSOLVER_PERMON): 			return_value = "PERMON QP solver"; break;
		case(GSOLVER_TAO): 				return_value = "TAO QP solver"; break;
	}
	return return_value;
}

template<class VectorBase>
std::string EntropyH1FEMModel<VectorBase>::print_thetasolvertype(ThetaSolverType thetasolvertype_in) const {
	std::string return_value = "?";
	switch(thetasolvertype_in){
		case(TSOLVER_AUTO): 			return_value = "AUTO"; break;
		case(TSOLVER_ENTROPY_DLIB): 	return_value = "ENTROPY_DLIB solver"; break;
		case(TSOLVER_ENTROPY_SPG): 		return_value = "ENTROPY_SPG solver"; break;
		case(TSOLVER_ENTROPY_NEWTON): 	return_value = "ENTROPY_NEWTON solver"; break;
	}
	return return_value;
}

/* constructor */
template<class VectorBase>
EntropyH1FEMModel<VectorBase>::EntropyH1FEMModel(TSData<VectorBase> &new_tsdata, int Km, double epssqr) {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyH1FEMModel<VectorBase>::~EntropyH1FEMModel(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T                 : " << this->tsdata->get_T() << std::endl;
	output <<  " - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output <<  " - scalef            : " << printbool(this->scalef) << std::endl;

	output <<  " - K                 : " << this->tsdata->get_K() << std::endl;
	output <<  " - Km                : " << this->get_Km() << std::endl;
	output <<  " - epssqr            : " << this->epssqr << std::endl;

	output <<  " - Graph             : " << std::endl;
	output.push();
	this->get_graph()->print(output);
	output.pop();

	output <<  " - thetalength       : " << this->thetavectorlength_global << std::endl;
	output <<  " - thetasolvertype   : " << print_thetasolvertype(this->thetasolvertype) << std::endl;
	output <<  " - gammasolvertype   : " << print_gammasolvertype(this->gammasolvertype) << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about model */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info" << std::endl;
	output_global <<  "  - T                 : " << this->tsdata->get_T() << std::endl;
	output_global <<  "  - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output_global <<  "  - scalef            : " << printbool(this->scalef) << std::endl;
	output_global <<  "  - K                 : " << this->tsdata->get_K() << std::endl;
	output_global <<  "  - Km                : " << this->get_Km() << std::endl;
	output_global <<  "  - epssqr            : " << this->epssqr << std::endl;
	output_global <<  "  - thetasolvertype   : " << print_thetasolvertype(this->thetasolvertype) << std::endl;
	output_global <<  "  - gammasolvertype   : " << print_gammasolvertype(this->gammasolvertype) << std::endl;

	output_global.push();
	this->get_graph()->print(output_global);
	output_global.pop();

	output_global <<  "  - thetalength       : " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables" << std::endl;
	output_global.push();
	output_local << "Tlocal =" << std::setw(6) << this->tsdata->get_Tlocal() << " (" << this->tsdata->get_Tbegin() << "," << this->tsdata->get_Tend() << "), ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

/* print model solution */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyH1FEMModel<VectorBase>::get_name() const {
	return "Entropy-H1-FEM Time-Series Model";
}

/* set new penalty */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::set_epssqr(double epssqr) {
	LOG_FUNC_BEGIN

	this->epssqr = epssqr;

	/* use old T to scale the function to obtain the same scale of function values (idea from Olga) */
	double coeff = this->epssqr;
	if(this->scalef){
		coeff *= (1.0/((double)(this->get_T())));
	} else {
		coeff *= 1.0;
	}

	if(this->A_shared != NULL){
		/* SPARSE */
		((BlockGraphSparseMatrix<VectorBase>*)A_shared)->set_coeff(coeff);
	}

	LOG_FUNC_END
}

/* prepare gamma solver */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::initialize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* prepare theta solver */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::initialize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destroy gamma solver */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

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
void EntropyH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

	LOG_FUNC_END
}

template<class VectorBase>
double EntropyH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver){
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

template<class VectorBase> //TODO: ?
void EntropyH1FEMModel<VectorBase>::get_linear_quadratic(double *linearL, double *quadraticL, GeneralSolver *gammasolver, GeneralSolver *thetasolver){
	*linearL = 0;
	*quadraticL = 0;
}

template<class VectorBase>
QPData<VectorBase>* EntropyH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
EntropyData<VectorBase>* EntropyH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}

template<class VectorBase>
BGMGraph<VectorBase> *EntropyH1FEMModel<VectorBase>::get_graph() const {
	return this->tsdata->get_decomposition()->get_graph();
}

template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::updateaftersolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


/* update theta solver */
template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void EntropyH1FEMModel<VectorBase>::updateaftersolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
double EntropyH1FEMModel<VectorBase>::get_aic(double L) const{
	return 2*log(L) + this->tsdata->get_K();
}

template<class VectorBase>
int EntropyH1FEMModel<VectorBase>::get_Km() const {
	return this->Km;
}

template<class VectorBase>
int EntropyH1FEMModel<VectorBase>::get_T() const {
	return this->tsdata->get_T();
}


}
} /* end namespace */

#endif
