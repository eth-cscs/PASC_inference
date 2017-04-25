#ifndef PASC_ENTROPYH1FEMMODEL_H
#define PASC_ENTROPYH1FEMMODEL_H

#ifndef USE_PETSC
 #error 'ENTROPYH1FEMMODEL is for PETSC'
#endif

#include "pascinference.h"

/* gamma problem */
#include "algebra/matrix/blockgraphsparse.h"

#include "data/qpdata.h"
#include "algebra/feasibleset/simplex_local.h"
#include "solver/spgqpsolver.h"
#include "solver/spgqpsolver_coeff.h"
#include "solver/taosolver.h"

#ifdef USE_PERMON
	#include "algebra/feasibleset/simplex_lineqbound.h"
	#include "solver/permonsolver.h"
#endif

/* theta problem */
#include "data/entropydata.h"
#include "solver/entropysolverdlib.h"
#include "solver/entropysolverspg.h"
#include "solver/entropysolvernewton.h"

#include "data/tsdata.h"


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
template<>
EntropyH1FEMModel<PetscVector>::EntropyH1FEMModel(TSData<PetscVector> &new_tsdata, int Km, double epssqr) {
	LOG_FUNC_BEGIN

	// TODO: enum in boost::program_options, not only int
	int gammasolvertype_int;
	int thetasolvertype_int;
	consoleArg.set_option_value("entropyh1femmodel_gammasolvertype", &gammasolvertype_int, GSOLVER_AUTO);
	consoleArg.set_option_value("entropyh1femmodel_thetasolvertype", &thetasolvertype_int, TSOLVER_AUTO);
	consoleArg.set_option_value("entropyh1femmodel_scalef", &scalef, GRAPHH1FEMMODEL_DEFAULT_SCALEF);
	
	this->gammasolvertype = static_cast<GammaSolverType>(gammasolvertype_int);
	this->thetasolvertype = static_cast<ThetaSolverType>(thetasolvertype_int);

	/* set given parameters */
	this->tsdata = &new_tsdata;
	this->Km = Km;
	
	/* prepare sequential vector with Theta - yes, all procesors will have the same information */
	this->thetavectorlength_local = tsdata->get_K()*this->get_Km();
	this->thetavectorlength_global = GlobalManager.get_size()*(this->thetavectorlength_local);

	/* set this model to data - tsdata will prepare gamma vector and thetavector */
	tsdata->set_model(*this);

	/* control the existence of graph */
	/* in this model, I need to work with graph, but maybe there is no graph in decomposition
	 * (for instance in the case without spatial decomposition), however it seems that user still wants to work 
	 * with this model based on graph. Therefore we create graph of disjoint nodes without any edge. This graph
	 * represents the situation and can be used for matrix-graph-based manipulation
	*/
	if(get_graph() == NULL){
		BGMGraph<PetscVector> *graph = this->tsdata->get_decomposition()->get_graph();
		double coordinates_array[this->tsdata->get_R()];

		for(int r=0;r<this->tsdata->get_R();r++){
			coordinates_array[r] = r;
		} 
		graph = new BGMGraph<PetscVector>(coordinates_array, this->tsdata->get_R(), 1);
		graph->process(0.0);
		
		this->tsdata->get_decomposition()->set_graph(*graph);
	}

	/* set regularization parameter */
	this->epssqr = epssqr;

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
template<>
void EntropyH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  "theta:" << std::endl;

	std::ostringstream temp;
	
	double *theta;
	TRYCXX( VecGetArray(thetadata->get_x()->get_vector(),&theta) );
	
	output_local.push();
	output_local << print_array(theta,tsdata->get_K(),get_Km()) << std::endl;
	output_local.pop();

	output_local.synchronize();
	output_global.synchronize();

	TRYCXX( VecRestoreArray(thetadata->get_x()->get_vector(),&theta) );

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
template<>
void EntropyH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* create data */
	gammadata = new QPData<PetscVector>();

	/* there is not reduction at all, we can use vectors from original data */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_x0(gammadata->get_x()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */

	/* moreover, for the residuum computation, we can use directly vector b */
	residuum = gammadata->get_b();

	/* use old T to scale the function to obtain the same scale of function values (idea from Olga) */
	double coeff = this->epssqr;
	if(scalef){
		coeff *= (1.0/((double)(this->get_T())));
	} else {
		coeff *= 1.0;
	}

	/* the vector of coefficient of blocks is set to NULL, therefore Theta will be not used to scale in penalisation */
	A_shared = new BlockGraphSparseMatrix<PetscVector>(*(this->tsdata->get_decomposition()), coeff, NULL );

	gammadata->set_A(A_shared); 

	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* automatic choice of solver */
	if(this->gammasolvertype == GSOLVER_AUTO){
		this->gammasolvertype = GSOLVER_SPGQP;
	}
	
	/* SPG-QP solver */
	if(this->gammasolvertype == GSOLVER_SPGQP){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolver<PetscVector>(*gammadata);
	}

	/* SPG-QP solver with special coefficient treatment */
	if(this->gammasolvertype == GSOLVER_SPGQP_COEFF){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolverC<PetscVector>(*gammadata);
	}

	/* Permon solver */
#ifdef USE_PERMON	
	if(this->gammasolvertype == GSOLVER_PERMON){
		/* the feasible set of QP is combination of linear equality constraints and bound inequality constraints */
		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(this->tsdata->get_decomposition()->get_T(),this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
		*gammasolver = new PermonSolver<PetscVector>(*gammadata);
	}
#endif

	/* TAO QP solver */
	if(this->gammasolvertype == GSOLVER_TAO){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(this->tsdata->get_decomposition()->get_T(),this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
		*gammasolver = new TaoSolver<PetscVector>(*gammadata);
	}

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	LOG_FUNC_END
}

/* prepare theta solver */
template<>
void EntropyH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN
	
	/* create data */
	thetadata = new EntropyData<PetscVector>(this->tsdata->get_T(), this->tsdata->get_K(), this->get_Km());
	thetadata->set_lambda(tsdata->get_thetavector());
	thetadata->set_x(tsdata->get_datavector());
	thetadata->set_gamma(tsdata->get_gammavector());
	thetadata->set_decomposition(tsdata->get_decomposition());

	/* automatic choice of solver */
	if(this->thetasolvertype == TSOLVER_AUTO){
		this->thetasolvertype = TSOLVER_ENTROPY_DLIB;
	}
	
	/* ENTROPY_DLIB solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_DLIB){
		/* create solver */
		*thetasolver = new EntropySolverDlib<PetscVector>(*thetadata);
	}

	/* ENTROPY_SPG solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_SPG){
		/* create solver */
		*thetasolver = new EntropySolverSPG<PetscVector>(*thetadata);
	}
	
	/* ENTROPY_NEWTON solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_NEWTON){
		/* create solver */
		*thetasolver = new EntropySolverNewton<PetscVector>(*thetadata);
	}	

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

template<>
void EntropyH1FEMModel<PetscVector>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	/* multiplicate vector b by coefficient */
	double coeff = -1.0;
	if(this->scalef){
		coeff *= (1.0/((double)(this->get_T())));
	}
	TRYCXX( VecScale(gammadata->get_b()->get_vector(), coeff) );

	LOG_FUNC_END
}

template<>
void EntropyH1FEMModel<PetscVector>::updateaftersolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* update theta solver */
template<>
void EntropyH1FEMModel<PetscVector>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<>
void EntropyH1FEMModel<PetscVector>::updateaftersolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	/* gammadata->get_b() = this->residuum */

	/* ENTROPY_DLIB solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_DLIB){
		/* create solver */
		((EntropySolverDlib<PetscVector> *)thetasolver)->compute_residuum(this->residuum); //TODO: retype?
	}

	/* ENTROPY_SPG solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_SPG){
		/* create solver */
		((EntropySolverSPG<PetscVector> *)thetasolver)->compute_residuum(this->residuum); //TODO: retype?
	}

	/* ENTROPY_NEWTON solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_NEWTON){
		/* create solver */
		((EntropySolverNewton<PetscVector> *)thetasolver)->compute_residuum(this->residuum); //TODO: retype?
	}

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
