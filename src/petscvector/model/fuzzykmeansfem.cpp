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

		/* create the residuum from original gamma vector */
		residuum = new GeneralVector<PetscVector>(*tsdata->get_gammavector()); //TODO: destroy this vector somewhere
		residuum_reduced = new GeneralVector<PetscVector>(*gammadata->get_x()); //TODO: destroy this vector somewhere

	} else {
		/* there is not reduction at all, we can use vectors from original data */
		gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */

		/* create the residuum from original gamma vector */
		residuum = new GeneralVector<PetscVector>(*tsdata->get_gammavector()); //TODO: destroy this vector somewhere
		residuum_reduced = residuum;
	}

	gamma_pow = new GeneralVector<PetscVector>(*tsdata->get_gammavector()); //TODO: destroy this vector somewhere

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

	int T = this->tsdata->get_decomposition()->get_T();
	int Tlocal = this->tsdata->get_decomposition()->get_Tlocal();

	int R = this->tsdata->get_decomposition()->get_R();
	int Rlocal = this->tsdata->get_decomposition()->get_Rlocal();

	int K = this->tsdata->get_decomposition()->get_K();

	int xdim = this->tsdata->get_decomposition()->get_xdim();

	/* if the problem is not reduced, then residuum=b, therefore it is not neccessary to perform reduction */
	if(fem->is_reduced()){
		this->fem->reduce_gamma(this->residuum, this->residuum_reduced);
		this->fem->reduce_gamma(tsdata->get_gammavector(), gammadata->get_x());
	}

    /* solve the problem on reduced mesh */
    double *residuum_reduced_arr, *gamma_reduced_arr;
	TRYCXX( VecGetArray(this->residuum_reduced->get_vector(), &residuum_reduced_arr) );
	TRYCXX( VecGetArray(gammadata->get_x()->get_vector(), &gamma_reduced_arr) );

	for(int t=0;t<fem->get_decomposition_reduced()->get_Tlocal();t++){
		for(int r=0;r<fem->get_decomposition_reduced()->get_Rlocal();r++){
			for(int k=0;k<K;k++){
                double sum = 0.0;
                double value;
                for(int k2=0;k2<K;k2++){
                    if(residuum_reduced_arr[t*K*Rlocal + r*K + k2] == 0){
                        value = 0;
                    } else {
                        value = residuum_reduced_arr[t*K*Rlocal + r*K + k]/residuum_reduced_arr[t*K*Rlocal + r*K + k2];
                    }
                    sum += pow(value, 1.0/((double)(fuzzifier - 1.0)) );
                }

                if(sum == 0){
                    gamma_reduced_arr[t*K*Rlocal + r*K + k] = 0.0;
                } else {
                    gamma_reduced_arr[t*K*Rlocal + r*K + k] = 1.0/sum;
                }
			}
		}
	}

	TRYCXX( VecRestoreArray(this->residuum_reduced->get_vector(), &residuum_reduced_arr) );
	TRYCXX( VecRestoreArray(gammadata->get_x()->get_vector(), &gamma_reduced_arr) );

	LOG_FUNC_END
}

template<>
void FuzzyKmeansModel<PetscVector>::gammasolver_updateaftersolve(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	/* if the problem is not reduced, then gammasolver->x = tsdata->gammavector, therefore it is not neccessary to perform prolongation */
	if(fem->is_reduced()){
		this->fem->prolongate_gamma(gammadata->get_x(), tsdata->get_gammavector());
	}

	LOG_FUNC_END
}


/* update theta solver */
template<>
void FuzzyKmeansModel<PetscVector>::thetasolver_updatebeforesolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	Vec gamma_Vec = tsdata->get_gammavector()->get_vector();
	Vec theta_Vec = tsdata->get_thetavector()->get_vector();
	Vec data_Vec = tsdata->get_datavector()->get_vector();

	/* now compute gamma^m */
	Vec gamma_pow_Vec = this->gamma_pow->get_vector();
	TRYCXX( VecCopy(gamma_Vec, gamma_pow_Vec) );
	TRYCXX( VecPow(gamma_pow_Vec, this->fuzzifier) );

	/* subvectors - data */
	Vec datan_Vec; /* vector corresponding to one dimension of data */
	IS datan_is;

	Vec gammak_Vec; /* gamma_k as a subvector of gamma */
	IS gammak_is;

	double gammakx;
	double gammaksum;

	/* get arrays */
	double *theta_arr;
	TRYCXX( VecGetArray(theta_Vec,&theta_arr) );

	int K = tsdata->get_K();
	int xdim = tsdata->get_xdim();

	double coeff = 1.0;

	/* through clusters */
	for(int k=0;k<K;k++){

		/* get gammak */
		this->tsdata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gamma_pow_Vec, gammak_is, &gammak_Vec) );

		/* compute gammaksum */
		TRYCXX( VecSum(gammak_Vec, &gammaksum) );

		for(int n=0;n<xdim;n++){

			/* get datan */
			this->tsdata->get_decomposition()->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_Vec, datan_is, &datan_Vec) );

			/* compute gammakx */
			TRYCXX( VecDot(datan_Vec, gammak_Vec, &gammakx) );

			if(gammaksum != 0){
				theta_arr[k*xdim + n] = gammakx/gammaksum;
			} else {
				theta_arr[k*xdim + n] = 0.0;
			}

			/* restore datan */
			TRYCXX( VecRestoreSubVector(data_Vec, datan_is, &datan_Vec) );
			TRYCXX( ISDestroy(&datan_is) );

		} /* endfor: through data dimension */

		TRYCXX( VecRestoreSubVector(gamma_pow_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	/* restore arrays */
	TRYCXX( VecRestoreArray(theta_Vec,&theta_arr) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

template<>
void FuzzyKmeansModel<PetscVector>::thetasolver_updateaftersolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

    /* compute new residuum for gamma problem */
    compute_residuum();

	LOG_FUNC_END
}

template<>
double FuzzyKmeansModel<PetscVector>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

    /* get residuum, here I suppose that this residuum has been already computed in theta step */
    Vec residuum_Vec = this->residuum->get_vector();

	/* now compute gamma^m */
	Vec gamma_Vec = tsdata->get_gammavector()->get_vector();
	Vec gamma_pow_Vec = this->gamma_pow->get_vector();
	TRYCXX( VecCopy(gamma_Vec, gamma_pow_Vec) );
	TRYCXX( VecPow(gamma_pow_Vec, this->fuzzifier) );

    double L;
    TRYCXX( VecDot(residuum_Vec, gamma_pow_Vec, &L) );

	LOG_FUNC_END

    return L;
}

template<>
void FuzzyKmeansModel<PetscVector>::compute_residuum(){
	LOG_FUNC_BEGIN

	int T = this->tsdata->get_decomposition()->get_T();
	int Tlocal = this->tsdata->get_decomposition()->get_Tlocal();

	int R = this->tsdata->get_decomposition()->get_R();
	int Rlocal = this->tsdata->get_decomposition()->get_Rlocal();

	int K = this->tsdata->get_decomposition()->get_K();

	int xdim = this->tsdata->get_decomposition()->get_xdim();

	/* compute residuum */
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

	/* restore arrays */
	TRYCXX( VecRestoreArray(residuum_Vec, &residuum_arr) );
	TRYCXX( VecRestoreArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	TRYCXX( VecRestoreArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );

	LOG_FUNC_END
}



}
} /* end namespace */

