#include "external/petscvector/model/entropyh1fem.h"

namespace pascinference {
namespace model {

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

	/* compute number of moments based on the xdim,Km, and K */

	coutMaster << "TEEEEEST: " << this->compute_number_of_moments() << std::endl;

	this->thetavectorlength_local = this->get_K()*this->compute_number_of_moments();

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
//		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(this->tsdata->get_decomposition()->get_T(),this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
//		*gammasolver = new PermonSolver<PetscVector>(*gammadata);
	}
#endif

	/* TAO QP solver */
	if(this->gammasolvertype == GSOLVER_TAO){
		/* the feasible set of QP is simplex */
//		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(this->tsdata->get_decomposition()->get_T(),this->tsdata->get_decomposition()->get_Tlocal(),this->tsdata->get_decomposition()->get_K())); 

		/* create solver */
//		*gammasolver = new TaoSolver<PetscVector>(*gammadata);
	}

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	LOG_FUNC_END
}

template<>
void EntropyH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	/* create data */
	thetadata = new EntropyData<PetscVector>(this->tsdata->get_T(), this->get_xdim(), this->tsdata->get_K(), this->get_Km());

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
	
	/* ENTROPY_NEWTON solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_NEWTON){
		/* create solver */
		*thetasolver = new EntropySolverNewton<PetscVector>(*thetadata);
	}	

	LOG_FUNC_END
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

template<>
void EntropyH1FEMModel<PetscVector>::updatebeforesolve_thetasolver(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	/* make sure that solved gamma is feasible */
	gammadata->get_feasibleset()->project(*(gammadata->get_x()));

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

	/* ENTROPY_NEWTON solver */
	if(this->thetasolvertype == TSOLVER_ENTROPY_NEWTON){
		/* create solver */
		((EntropySolverNewton<PetscVector> *)thetasolver)->compute_residuum(this->residuum); //TODO: retype?
	}

	LOG_FUNC_END
}



}
} /* end namespace */

