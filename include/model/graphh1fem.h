#ifndef PASC_GRAPHH1FEMMODEL_H
#define PASC_GRAPHH1FEMMODEL_H

#ifndef USE_PETSCVECTOR
 #error 'GRAPHH1FEMMODEL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

/* for reduction and prolongation of data */
#include "common/fem.h"
#include "common/femhat.h"

/* gamma problem */
//#include "algebra/matrix/blockgraphfree.h" // TODO: implement?
#include "algebra/matrix/blockgraphsparse.h"

#include "algebra/feasibleset/simplex_local.h"
#include "solver/spgqpsolver.h"
#include "solver/spgqpsolver_coeff.h"
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
			SOLVER_PERMON=3 /**< PERMONQP solver (augumented lagrangians combined with active-set method) */
		} GammaSolverType;

		/** @brief type of FEM used to reduce gamma problem 
		 * 
		 */
		typedef enum { 
			FEM_SUM=0, /**< use the sum */
			FEM_HAT=1 /**< use linear "hat" functions */
		} FemType;

	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex, will be solved by SPG-QP */
	 	SimpleData<VectorBase> *thetadata; /**< this problem is solved during assembly  */

		double epssqr; /**< penalty coeficient */
		
		/* for theta problem */
		GeneralMatrix<VectorBase> *A_shared; /**< matrix shared by gamma and theta solver */
		GeneralVector<VectorBase> *Agamma; /**< temp vector for storing A_shared*gamma */
		GeneralVector<VectorBase> *residuum; /**< temp vector for residuum computation */
		
		bool usethetainpenalty; /**< use the value of Theta in penalty parameter to scale blocks */
		
		GammaSolverType gammasolvertype; /**< the type of used solver */
		
		double fem_reduce; /**< coeficient of mesh reduction */
		int T_reduced; /**< the length of reduced gamma problem, T_reduced = fem_double*T */
		Decomposition *decomposition_reduced; /**< the decomposition of reduced problem */
		Fem *fem;
		FemType femtype; /**< the type of used FEM for reduction */

	public:

		/** @brief constructor from data and regularisation constant
		 * 
		 * @param tsdata time-series data on which model operates
		 * @param epssqr regularisation constant
		 */ 	
		GraphH1FEMModel(TSData<VectorBase> &tsdata, double epssqr, double fem_reduce = 1.0, bool usethetainpenalty = false);

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
		BGMGraph *get_graph() const;

		GeneralVector<VectorBase> *get_coordinatesVTK() const;
		int get_coordinatesVTK_dim() const;

		double get_aic(double L) const;
		
		void set_epssqr(double epssqr);
		bool get_usethetainpenalty() const;
		
		int get_T_reduced() const;
		int get_T() const;
		Decomposition *get_decomposition_reduced() const;
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
GraphH1FEMModel<PetscVector>::GraphH1FEMModel(TSData<PetscVector> &new_tsdata, double epssqr, double fem_reduce, bool usethetainpenalty) {
	LOG_FUNC_BEGIN

	// TODO: enum in boost::program_options, not only int
	int gammasolvertype_int;
	consoleArg.set_option_value("graphh1femmodel_gammasolvertype", &gammasolvertype_int, SOLVER_AUTO);
	this->gammasolvertype = static_cast<GammaSolverType>(gammasolvertype_int);

	// TODO: enum in boost::program_options, not only int
	int femtype_int;
	consoleArg.set_option_value("graphh1femmodel_femtype", &femtype_int, FEM_HAT);
	this->femtype = static_cast<FemType>(femtype_int);

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
		BGMGraph *graph = this->tsdata->get_decomposition()->get_graph();
		double coordinates_array[this->tsdata->get_R()];

		for(int r=0;r<this->tsdata->get_R();r++){
			coordinates_array[r] = r;
		} 
		graph = new BGMGraph(coordinates_array, this->tsdata->get_R(), 1);
		graph->process(0.0);
		
		this->tsdata->get_decomposition()->set_graph(*graph);
	}

	/* prepare parameters of reduced problem */
	this->fem_reduce = fem_reduce;
	this->T_reduced = ceil(this->tsdata->get_T()*fem_reduce);
	if(fem_reduce < 1.0){
		/* compute new decomposition */
		this->decomposition_reduced = new Decomposition(this->T_reduced, 
				*(this->tsdata->get_decomposition()->get_graph()), 
				this->tsdata->get_decomposition()->get_K(), 
				this->tsdata->get_decomposition()->get_xdim(), 
				this->tsdata->get_decomposition()->get_DDT_size(), 
				this->tsdata->get_decomposition()->get_DDR_size());
		
		/* set the type of FEM */
		if(this->femtype == FEM_SUM){
			fem = new Fem(this->tsdata->get_decomposition(),this->decomposition_reduced);
		}
		if(this->femtype == FEM_HAT){
			fem = new FemHat(this->tsdata->get_decomposition(),this->decomposition_reduced);
		}

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->decomposition_reduced = this->tsdata->get_decomposition();
	}

	/* set regularization parameter */
	this->set_epssqr(epssqr);
	
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

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T                 : " << this->tsdata->get_T() << std::endl;
	output <<  " - xdim              : " << this->tsdata->get_xdim() << std::endl;

	/* information of reduced problem */
	output <<  " - fem_reduce        : " << this->fem_reduce << std::endl;
	output <<  "  - femtype          : ";
	switch(this->femtype){
		case(FEM_SUM): output << "FEM_SUM"; break;
		case(SOLVER_SPGQP): output << "FEM_HAT"; break;
	}
	output << std::endl;
	output <<  "  - T_reduced        : " << this->T_reduced << std::endl;
	output <<  "  - reduced decomposition : " << std::endl;
	output.push();
	output.push();
	decomposition_reduced->print(output);
	output.pop();
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

	/* information of reduced problem */
	output_global <<  "  - fem_reduce        : " << this->fem_reduce << std::endl;
	output_global <<  "   - femtype          : ";
	switch(this->femtype){
		case(FEM_SUM): output_global << "FEM_SUM"; break;
		case(SOLVER_SPGQP): output_global << "FEM_HAT"; break;
	}
	output_global << std::endl;
	output_global <<  "   - T_reduced        : " << this->T_reduced << std::endl;
	output_global <<  "   - reduced decomposition : " << std::endl;
	output_global.push();
	output_global.push();
	decomposition_reduced->print_content(output_global,output_local);
	output_global.pop();
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
	this->epssqr = epssqr;

	/* use old T to scale the function to obtain the same scale of function values (idea from Olga) */
//	double coeff = (1.0/((double)(this->get_T_reduced())))*this->epssqr;
	double coeff = (1.0/((double)(this->get_T())))*this->epssqr;

	if(this->A_shared){
		/* SPARSE */
		((BlockGraphSparseMatrix<VectorBase>*)A_shared)->set_coeff(coeff);
	}
}

/* prepare gamma solver */
template<>
void GraphH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* create data */
	gammadata = new QPData<PetscVector>();

	/* deal with problem reduction */
	if(this->fem_reduce < 1.0){
		/* there is a reduction, we have to create new reduced gammavector */
		Vec x_reduced_Vec;
		this->decomposition_reduced->createGlobalVec_gamma(&x_reduced_Vec);
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
	double coeff = (1.0/((double)(this->get_T())))*this->epssqr;

	/* SPARSE */
	if(usethetainpenalty){
		/* use thetavector as a vector of coefficient for scaling blocks */
		A_shared = new BlockGraphSparseMatrix<PetscVector>(*(this->decomposition_reduced), coeff, tsdata->get_thetavector() );
	} else {
		/* the vector of coefficient of blocks is set to NULL, therefore Theta will be not used to scale in penalisation */
		A_shared = new BlockGraphSparseMatrix<PetscVector>(*(this->decomposition_reduced), coeff, NULL );
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
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(this->decomposition_reduced->get_Tlocal()*this->decomposition_reduced->get_Rlocal(),this->decomposition_reduced->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolver<PetscVector>(*gammadata);

		/* modify stopping criteria based on reduction */
		if(this->fem_reduce < 1.0){
//			(*gammasolver)->set_eps((this->fem_reduce)*(*gammasolver)->get_eps());
		}
	}

	/* SPG-QP solver */
	if(this->gammasolvertype == SOLVER_SPGQP_COEFF){
		/* the feasible set of QP is simplex */
		gammadata->set_feasibleset(new SimplexFeasibleSet_Local<PetscVector>(this->decomposition_reduced->get_Tlocal()*this->decomposition_reduced->get_Rlocal(),this->decomposition_reduced->get_K())); 

		/* create solver */
		*gammasolver = new SPGQPSolverC<PetscVector>(*gammadata);

		/* modify stopping criteria based on reduction */
		if(this->fem_reduce < 1.0){
//			(*gammasolver)->set_eps((this->fem_reduce)*(*gammasolver)->get_eps());
		}
	}

	/* Permon solver */
#ifdef USE_PERMON	
	if(this->gammasolvertype == SOLVER_PERMON){
		/* the feasible set of QP is combination of linear equality constraints and bound inequality constraints */
		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(this->decomposition_reduced->get_T()*this->decomposition_reduced->get_R(),this->decomposition_reduced->get_Tlocal()*this->decomposition_reduced->get_Rlocal(),this->decomposition_reduced->get_K())); 

		/* create solver */
		*gammasolver = new PermonSolver<PetscVector>(*gammadata);

		/* modify stopping criteria based on reduction */
		if(this->fem_reduce < 1.0){
//			(*gammasolver)->set_eps((this->fem_reduce)*(*gammasolver)->get_eps());
		}
	}
#endif

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
	if(this->fem_reduce < 1){
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
BGMGraph *GraphH1FEMModel<VectorBase>::get_graph() const {
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
	if(fem_reduce < 1.0){
		this->fem->reduce_gamma(this->residuum, gammadata->get_b());
		this->fem->reduce_gamma(tsdata->get_gammavector(), gammadata->get_x());
	}

	/* multiplicate vector b by coefficient */
	double coeff = (-1.0/((double)(this->get_T_reduced())));
//	double coeff = (-1.0/((double)(this->get_T())));
	TRYCXX( VecScale(gammadata->get_b()->get_vector(), coeff) );

	LOG_FUNC_END
}

template<>
void GraphH1FEMModel<PetscVector>::updateaftersolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	/* if the problem is not reduced, then gammasolver->x = tsdata->gammavector, therefore it is not neccessary to perform prolongation */
	if(fem_reduce < 1.0){
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
	double return_value;
	if(this->fem_reduce < 1.0){
		return_value = this->T_reduced;
	} else {
		return_value = this->get_T();
	}

	return return_value;
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
Decomposition *GraphH1FEMModel<VectorBase>::get_decomposition_reduced() const {
	return this->decomposition_reduced;
}



}
} /* end namespace */

#endif
