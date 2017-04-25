#ifndef PASC_GRAPHH1FEMMODEL_H
#define PASC_GRAPHH1FEMMODEL_H

#ifndef USE_PETSC
 #error 'GRAPHH1FEMMODEL is for PETSC'
#endif

#include "pascinference.h"

/* for reduction and prolongation of data */
#include "common/fem.h"

/* gamma problem */
//#include "algebra/matrix/blockgraphfree.h" // TODO: implement?
#include "algebra/matrix/blockgraphsparse.h"

#include "algebra/feasibleset/simplex_local.h"
#include "solver/spgqpsolver.h"
#include "solver/spgqpsolver_coeff.h"
#include "solver/taosolver.h"
#include "data/qpdata.h"

#ifdef USE_PERMON
	#include "solver/permonsolver.h"
	#include "algebra/feasibleset/simplex_lineqbound.h"
#endif

/* theta problem */
#include "solver/simplesolver.h"
#include "data/simpledata.h"

#include "data/tsdata.h"


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
template<>
GraphH1FEMModel<PetscVector>::GraphH1FEMModel(TSData<PetscVector> &new_tsdata, double epssqr, Fem<PetscVector> *new_fem, bool usethetainpenalty) {
	LOG_FUNC_BEGIN

	// TODO: enum in boost::program_options, not only int
	int gammasolvertype_int;
	consoleArg.set_option_value("graphh1femmodel_gammasolvertype", &gammasolvertype_int, SOLVER_AUTO);
	consoleArg.set_option_value("graphh1femmodel_scalef", &scalef, GRAPHH1FEMMODEL_DEFAULT_SCALEF);
	
	this->gammasolvertype = static_cast<GammaSolverType>(gammasolvertype_int);

	/* set given parameters */
	this->usethetainpenalty = usethetainpenalty;
	this->tsdata = &new_tsdata;
	
	/* prepare sequential vector with Theta - yes, all procesors will have the same information */
	this->thetavectorlength_local = tsdata->get_K()*tsdata->get_xdim();
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

	/* prepare parameters and decomposition of reduced problem */
	/* if FEM is not given, then prepare FEM without reduction */
	if(new_fem == NULL){
		this->fem = new Fem<PetscVector>(1.0);
	} else {
		this->fem = new_fem;
	}
	
	double fem_reduce = this->fem->get_fem_reduce();

	this->fem->set_decomposition_original(this->tsdata->get_decomposition());
	this->fem->compute_decomposition_reduced();

	/* set regularization parameter */
	this->epssqr = epssqr;

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
template<>
void GraphH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  "theta:" << std::endl;

	std::ostringstream temp;
	
	double *theta;
	TRYCXX( VecGetArray(thetadata->get_x()->get_vector(),&theta) );
	
	int k,n;

	output_local.push();
	output_local << "- proc: " << GlobalManager.get_rank() << std::endl;
	output_local.push();
	for(k=0;k<tsdata->get_K();k++){
		output_local <<  "- k = " << k << std::endl;

		/* mu */
		output_local.push();
		output_local <<  "- theta = [";
		for(n=0;n<tsdata->get_xdim();n++){
			temp << theta[k*tsdata->get_xdim() + n];
			output_local << temp.str();
			if(n < tsdata->get_xdim()-1){
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
template<>
void GraphH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* create data */
	gammadata = new QPData<PetscVector>();

	/* deal with problem reduction */
	if(fem->is_reduced()){
		/* there is a reduction, we have to create new reduced gammavector */
		Vec x_reduced_Vec;
		get_decomposition_reduced()->createGlobalVec_gamma(&x_reduced_Vec);
		gammadata->set_x(new GeneralVector<PetscVector>(x_reduced_Vec)); /* create new linear term of QP problem */
		gammadata->set_x0(gammadata->get_x()); /* the initial approximation of QP problem is gammavector */
		gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */
		
		/* create the residuum from original gamma vector */
		residuum = new GeneralVector<PetscVector>(*tsdata->get_gammavector());
	
	} else {
		/* there is not reduction at all, we can use vectors from original data */
		gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
		gammadata->set_x0(gammadata->get_x()); /* the initial approximation of QP problem is gammavector */
		gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */

		/* moreover, for the residuum computation, we can use directly vector b */
		residuum = gammadata->get_b();
	}

	/* use old T to scale the function to obtain the same scale of function values (idea from Olga) */
//	double coeff = (1.0/((double)(this->get_T_reduced())))*this->epssqr;
	double coeff = this->epssqr;
	if(scalef){
		coeff *= (1.0/((double)(this->get_T())));
	} else {
		coeff *= ((double)(this->get_T_reduced())/((double)(this->get_T())));
	}

	/* SPARSE */
	if(usethetainpenalty){
		/* use thetavector as a vector of coefficient for scaling blocks */
		A_shared = new BlockGraphSparseMatrix<PetscVector>(*(get_decomposition_reduced()), coeff, tsdata->get_thetavector() );
	} else {
		/* the vector of coefficient of blocks is set to NULL, therefore Theta will be not used to scale in penalisation */
		A_shared = new BlockGraphSparseMatrix<PetscVector>(*(get_decomposition_reduced()), coeff, NULL );
	}

	gammadata->set_A(A_shared); 

	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
//	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	/* automatic choice of solver */
	if(this->gammasolvertype == SOLVER_AUTO){
		this->gammasolvertype = SOLVER_SPGQP;
	}
	
	/* SPG-QP solver */
	if(this->gammasolvertype == SOLVER_SPGQP){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(get_decomposition_reduced()->get_Tlocal()*get_decomposition_reduced()->get_Rlocal(),get_decomposition_reduced()->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolver<PetscVector>(*gammadata);
	}

	/* SPG-QP solver with special coefficient treatment */
	if(this->gammasolvertype == SOLVER_SPGQP_COEFF){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(get_decomposition_reduced()->get_Tlocal()*get_decomposition_reduced()->get_Rlocal(),get_decomposition_reduced()->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolverC<PetscVector>(*gammadata);
	}

	/* Permon solver */
#ifdef USE_PERMON	
	if(this->gammasolvertype == SOLVER_PERMON){
		/* the feasible set of QP is combination of linear equality constraints and bound inequality constraints */
		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(get_decomposition_reduced()->get_T()*get_decomposition_reduced()->get_R(),get_decomposition_reduced()->get_Tlocal()*get_decomposition_reduced()->get_Rlocal(),get_decomposition_reduced()->get_K())); 

		/* create solver */
		*gammasolver = new PermonSolver<PetscVector>(*gammadata);
	}
#endif

	/* TAO QP solver */
	if(this->gammasolvertype == SOLVER_TAO){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(get_decomposition_reduced()->get_T()*get_decomposition_reduced()->get_R(),get_decomposition_reduced()->get_Tlocal()*get_decomposition_reduced()->get_Rlocal(),get_decomposition_reduced()->get_K())); 

		/* create solver */
		*gammasolver = new TaoSolver<PetscVector>(*gammadata);
	}

	LOG_FUNC_END
}

/* prepare theta solver */
template<>
void GraphH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN
	
	/* create data */
	thetadata = new SimpleData<PetscVector>();
	thetadata->set_x(tsdata->get_thetavector());

	/* create solver */
	*thetasolver = new SimpleSolver<PetscVector>(*thetadata);
	
	/* create aux vector for gamma^T A gamma */
	Vec Agamma_Vec;
	TRYCXX( VecDuplicate(tsdata->get_gammavector()->get_vector(),&Agamma_Vec) );
	Agamma = new GeneralVector<PetscVector>(Agamma_Vec);
	
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

template<>
void GraphH1FEMModel<PetscVector>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	int T = this->tsdata->get_decomposition()->get_T();
	int Tlocal = this->tsdata->get_decomposition()->get_Tlocal();

	int R = this->tsdata->get_decomposition()->get_R();
	int Rlocal = this->tsdata->get_decomposition()->get_Rlocal();

	int K = this->tsdata->get_decomposition()->get_K();

	/* update gamma_solver data - prepare new linear term */
	const double *theta_arr;
	TRYCXX( VecGetArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	
	const double *data_arr;
	TRYCXX( VecGetArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	
	double *residuum_arr;
	TRYCXX( VecGetArray(this->residuum->get_vector(), &residuum_arr) );

	for(int t=0;t<Tlocal;t++){
		for(int r=0;r<Rlocal;r++){
			for(int k=0;k<K;k++){
				residuum_arr[t*K*Rlocal + r*K + k] = (data_arr[t*Rlocal+r] - theta_arr[k])*(data_arr[t*Rlocal+r] - theta_arr[k]);
			}
		}
	}

	/* coeffs of A_shared are updated via computation of Theta :) */

	/* restore arrays */
	TRYCXX( VecRestoreArray(this->residuum->get_vector(), &residuum_arr) );
	TRYCXX( VecRestoreArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	TRYCXX( VecRestoreArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );

	/* if the problem is not reduced, then residuum=b, therefore it is not neccessary to perform reduction */
	if(fem->is_reduced()){
		this->fem->reduce_gamma(this->residuum, gammadata->get_b());
		this->fem->reduce_gamma(tsdata->get_gammavector(), gammadata->get_x());
	}

	/* multiplicate vector b by coefficient */
	double coeff = -1.0;
	if(this->scalef){
		coeff *= (1.0/((double)(this->get_T_reduced())));
	}
	
//	double coeff = (-1.0/((double)(this->get_T())));
	TRYCXX( VecScale(gammadata->get_b()->get_vector(), coeff) );

	LOG_FUNC_END
}

template<>
void GraphH1FEMModel<PetscVector>::updateaftersolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	/* if the problem is not reduced, then gammasolver->x = tsdata->gammavector, therefore it is not neccessary to perform prolongation */
	if(fem->is_reduced()){
		this->fem->prolongate_gamma(gammadata->get_x(), tsdata->get_gammavector());
	}

	LOG_FUNC_END
}


/* update theta solver */
template<>
void GraphH1FEMModel<PetscVector>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	// TODO: if Theta is not in penalty term, this computation is completely WRONG! However, it doesn't matter if Theta is given

	Vec gamma_Vec = tsdata->get_gammavector()->get_vector();
	Vec theta_Vec = tsdata->get_thetavector()->get_vector();
	Vec data_Vec = tsdata->get_datavector()->get_vector();

	/* I will use A_shared with coefficients equal to 1, therefore I set Theta=1 */
	TRYCXX( VecSet(tsdata->get_thetavector()->get_vector(),1.0) );
	TRYCXX( VecAssemblyBegin(tsdata->get_thetavector()->get_vector()) );
	TRYCXX( VecAssemblyEnd(tsdata->get_thetavector()->get_vector()) );

	/* now compute A*gamma */
	Vec Agamma_Vec;
	if(usethetainpenalty){
		/* only if Theta is in penalty term */
		*Agamma = (*A_shared)*(*(tsdata->get_gammavector()));
		Agamma_Vec = Agamma->get_vector();
	}

	/* subvectors */
	Vec gammak_Vec;
	Vec Agammak_Vec;
	IS gammak_is;
	
	double gammakAgammak;
	double gammakx;
	double gammaksum;
	
	/* get arrays */
	double *theta_arr;
	TRYCXX( VecGetArray(theta_Vec,&theta_arr) );

	int K = tsdata->get_K();

	double coeff = 1.0;
//	double coeff = 1.0/((double)(tsdata->get_R()*tsdata->get_T()));

	/* through clusters */
	for(int k=0;k<K;k++){
		
		/* get gammak */
		this->tsdata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		/* compute gammakAgammak */
		if(usethetainpenalty){
			/* only if Theta is in penalty term */
			TRYCXX( VecGetSubVector(Agamma_Vec, gammak_is, &Agammak_Vec) );
			TRYCXX( VecDot(gammak_Vec, Agammak_Vec, &gammakAgammak) );
		}
		
		/* compute gammakx */
		TRYCXX( VecDot(data_Vec, gammak_Vec, &gammakx) );

		/* compute gammaksum */
		TRYCXX( VecSum(gammak_Vec, &gammaksum) );

		if(usethetainpenalty){
			/* only if Theta is in penalty term */
			if(coeff*gammaksum + 0.5*gammakAgammak != 0){
				theta_arr[k] = (coeff*gammakx)/(coeff*gammaksum + 0.5*gammakAgammak);
			} else {
				theta_arr[k] = 0.0;
			}
		} else {
			/* if Theta is not in penalty term, then the computation is based on kmeans */
			if(gammaksum != 0){
				theta_arr[k] = gammakx/gammaksum;
			} else {
				theta_arr[k] = 0.0;
			}
		}
	
		TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		if(usethetainpenalty){
			/* only if Theta is in penalty term */
			TRYCXX( VecRestoreSubVector(Agamma_Vec, gammak_is, &Agammak_Vec) );
		}
		TRYCXX( ISDestroy(&gammak_is) );
	}	

	/* restore arrays */
	TRYCXX( VecRestoreArray(theta_Vec,&theta_arr) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

template<>
void GraphH1FEMModel<PetscVector>::updateaftersolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

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
