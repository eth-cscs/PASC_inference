#include "external/petscvector/model/graphh1fem.h"

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

	output_local.push();
	output_local << "- proc: " << GlobalManager.get_rank() << std::endl;
	output_local.push();
	for(int k=0;k<tsdata->get_K();k++){
		output_local <<  "- k = " << k << std::endl;

		/* mu */
		output_local.push();
		output_local <<  "- theta = [";
		for(int n=0;n<tsdata->get_xdim();n++){
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
//		gammadata->set_feasibleset(new SimplexFeasibleSet_LinEqBound<PetscVector>(get_decomposition_reduced()->get_T()*get_decomposition_reduced()->get_R(),get_decomposition_reduced()->get_Tlocal()*get_decomposition_reduced()->get_Rlocal(),get_decomposition_reduced()->get_K()));

		/* create solver */
//		*gammasolver = new TaoSolver<PetscVector>(*gammadata);
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

template<>
void GraphH1FEMModel<PetscVector>::updatebeforesolve_gammasolver(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	int T = this->tsdata->get_decomposition()->get_T();
	int Tlocal = this->tsdata->get_decomposition()->get_Tlocal();

	int R = this->tsdata->get_decomposition()->get_R();
	int Rlocal = this->tsdata->get_decomposition()->get_Rlocal();

	int K = this->tsdata->get_decomposition()->get_K();

	int xdim = this->tsdata->get_decomposition()->get_xdim();


	/* update gamma_solver data - prepare new linear term */
	const double *theta_arr;
	TRYCXX( VecGetArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );

	const double *data_arr;
	TRYCXX( VecGetArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );

	double *residuum_arr;
	Vec residuum_Vec = this->residuum->get_vector();
	TRYCXX( VecSet(residuum_Vec,0.0) );
	TRYCXX( VecGetArray(residuum_Vec, &residuum_arr) );

	for(int t=0;t<Tlocal;t++){
		for(int r=0;r<Rlocal;r++){
			for(int k=0;k<K;k++){
				for(int n=0;n<xdim;n++){
					residuum_arr[t*K*Rlocal + r*K + k] += (data_arr[(t*Rlocal+r)*xdim + n] - theta_arr[k*xdim+n])*(data_arr[(t*Rlocal+r)*xdim + n] - theta_arr[k*xdim+n]);
				}
			}
		}
	}

	/* coeffs of A_shared are updated via computation of Theta :) */

	/* restore arrays */
	TRYCXX( VecRestoreArray(residuum_Vec, &residuum_arr) );
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

	/* subvectors - data */
	Vec datan_Vec; /* vector corresponding to one dimension of data */
	IS datan_is;

	Vec gammak_Vec; /* gamma_k as a subvector of gamma */
	Vec Agammak_Vec;
	IS gammak_is;

	double gammakAgammak;
	double gammakx;
	double gammaksum;

	/* get arrays */
	double *theta_arr;
	TRYCXX( VecGetArray(theta_Vec,&theta_arr) );

	int K = tsdata->get_K();
	int xdim = tsdata->get_xdim();

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

		/* compute gammaksum */
		TRYCXX( VecSum(gammak_Vec, &gammaksum) );

		for(int n=0;n<xdim;n++){

			/* get datan */
			this->tsdata->get_decomposition()->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_Vec, datan_is, &datan_Vec) );

			/* compute gammakx */
			TRYCXX( VecDot(datan_Vec, gammak_Vec, &gammakx) );

			if(usethetainpenalty){
				/* only if Theta is in penalty term */
				if(coeff*gammaksum + 0.5*gammakAgammak != 0){
					theta_arr[k*xdim + n] = (coeff*gammakx)/(coeff*gammaksum + 0.5*gammakAgammak);
				} else {
					theta_arr[k*xdim + n] = 0.0;
				}
			} else {
				/* if Theta is not in penalty term, then the computation is based on kmeans */
				if(gammaksum != 0){
					theta_arr[k*xdim + n] = gammakx/gammaksum;
				} else {
					theta_arr[k*xdim + n] = 0.0;
				}
			}

			/* restore datan */
			TRYCXX( VecRestoreSubVector(data_Vec, datan_is, &datan_Vec) );
			TRYCXX( ISDestroy(&datan_is) );

		} /* endfor: through data dimension */

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

}
} /* end namespace */

