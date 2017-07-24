#include "external/petscvector/model/fuzzykmeans.h"

namespace pascinference {
namespace model {

/* constructor */
template<>
FuzzyKmeansModel<PetscVector>::FuzzyKmeansModel(TSData<PetscVector> &new_tsdata, double new_fuzzifier, Fem<PetscVector> *new_fem) {
	LOG_FUNC_BEGIN

	/* set given parameters */
	this->tsdata = &new_tsdata;
	this->fuzzifier = new_fuzzifier;

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

	LOG_FUNC_END
}

/* print model solution */
template<>
void FuzzyKmeansModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
void FuzzyKmeansModel<PetscVector>::gammasolver_initialize(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* create data */
	gammadata = new SimpleData<PetscVector>();

	/* deal with problem reduction */
	if(fem->is_reduced()){
		/* there is a reduction, we have to create new reduced gammavector */
		Vec x_reduced_Vec;
		get_decomposition_reduced()->createGlobalVec_gamma(&x_reduced_Vec);
		gammadata->set_x(new GeneralVector<PetscVector>(x_reduced_Vec)); /* create new linear term of QP problem */
	} else {
		/* there is not reduction at all, we can use vectors from original data */
		gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	}

	/* create the residuum from original gamma vector */
	residuum = new GeneralVector<PetscVector>(*tsdata->get_gammavector());

	/* create solver */
	*gammasolver = new SimpleSolver<PetscVector>(*gammadata);

	LOG_FUNC_END
}

/* prepare theta solver */
template<>
void FuzzyKmeansModel<PetscVector>::thetasolver_initialize(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	/* create data */
	thetadata = new SimpleData<PetscVector>();
	thetadata->set_x(tsdata->get_thetavector());

	/* create solver */
	*thetasolver = new SimpleSolver<PetscVector>(*thetadata);

	LOG_FUNC_END
}

template<>
void FuzzyKmeansModel<PetscVector>::gammasolver_updatebeforesolve(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<>
void FuzzyKmeansModel<PetscVector>::gammasolver_updateaftersolve(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* update theta solver */
template<>
void FuzzyKmeansModel<PetscVector>::thetasolver_updatebeforesolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<>
void FuzzyKmeansModel<PetscVector>::thetasolver_updateaftersolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

}
} /* end namespace */

