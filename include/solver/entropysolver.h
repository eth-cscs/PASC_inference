/** @file entropysolver.h
 *  @brief Solver which solves problem with integrals from Anna
 *
 *  @author Anna Marchenko & Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVER_H
#define	PASC_ENTROPYSOLVER_H

#include "pascinference.h"
#include "data/entropydata.h"

namespace pascinference {
namespace solver {

/** \class EntropySolver
 *  \brief Solver which solves problem with integrals from Anna
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolver: public GeneralSolver {
	protected:
		Timer timer_solve; /**< total solution time */
		Timer timer_compute_moments; /**< time of moment computation */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments */
		GeneralVector<VectorBase> *x_power; /**< temp vector for storing power of x */
		GeneralVector<VectorBase> *x_power_gammak; /**< temp vector for storing power of x * gamma_k */

	public:

		EntropySolver();
		EntropySolver(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolver();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printcontent(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		EntropyData<VectorBase> *get_data() const;

		void compute_moments();
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
EntropySolver<VectorBase>::EntropySolver(){
	LOG_FUNC_BEGIN
	
	entropydata = NULL;
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	LOG_FUNC_END
}

template<class VectorBase>
EntropySolver<VectorBase>::EntropySolver(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* prepare timers */
	this->timer_solve.restart();	

	/* prepare auxiliary vectors */
	x_power = new GeneralVector<PetscVector>(*entropydata->get_x());
	x_power_gammak = new GeneralVector<PetscVector>(*entropydata->get_x());
	
	/* create aux vector for the computation of moments */
	Vec moments_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&moments_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(moments_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(moments_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(moments_Vec,entropydata->get_K()*entropydata->get_Km(),PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(moments_Vec) );
	this->moments = new GeneralVector<PetscVector>(moments_Vec);
	

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolver<VectorBase>::~EntropySolver(){
	LOG_FUNC_BEGIN

	free(x_power);
	free(x_power_gammak);
	free(moments);

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		coutMaster.push();
		entropydata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void EntropySolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << std::endl;
	output_global <<  " - eps:        " << this->eps << std::endl;
	output_global <<  " - debugmode: " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		output_global << "- data:" << std::endl;
		coutMaster.push();
		entropydata->print(output_global,output_local);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void EntropySolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(entropydata){
		output << "- data:" << std::endl;
		coutMaster.push();
		entropydata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve   = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_moments = " << this->timer_compute_moments.get_value_sum() << std::endl;


	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolver<VectorBase>::get_name() const {
	std::string return_value = "EntropySolver<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
EntropyData<VectorBase> *EntropySolver<VectorBase>::get_data() const {
	return entropydata;
}

template<class VectorBase>
void EntropySolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	this->compute_moments();
	
	coutMaster << "Moments: " << *moments << std::endl;

	this->timer_solve.start(); 


	///* this solve does not anything */
	//Vec gamma_Vec = tsdata->get_gammavector()->get_vector();
	//Vec theta_Vec = tsdata->get_thetavector()->get_vector();
	//Vec data_Vec = tsdata->get_datavector()->get_vector();

	///* I will use A_shared with coefficients equal to 1, therefore I set Theta=1 */
	//TRYCXX( VecSet(tsdata->get_thetavector()->get_vector(),1.0) );
	//TRYCXX( VecAssemblyBegin(tsdata->get_thetavector()->get_vector()) );
	//TRYCXX( VecAssemblyEnd(tsdata->get_thetavector()->get_vector()) );

	///* subvectors */
	//Vec gammak_Vec;
	//Vec Agammak_Vec;
	//IS gammak_is;
	
	//double gammakx;
	//double gammaksum;
	
	///* get arrays */
	//double *theta_arr;
	//TRYCXX( VecGetArray(theta_Vec,&theta_arr) );

	//int K = tsdata->get_K();

	//double coeff = 1.0;

	///* through clusters */
	//for(int k=0;k<K;k++){
		
		///* get gammak */
		//this->tsdata->get_decomposition()->createIS_gammaK(&gammak_is, k);
		//TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		///* compute gammakx */
		//TRYCXX( VecDot(data_Vec, gammak_Vec, &gammakx) );

		///* compute gammaksum */
		//TRYCXX( VecSum(gammak_Vec, &gammaksum) );

		///* if Theta is not in penalty term, then the computation is based on kmeans */
		//if(gammaksum != 0){
			//theta_arr[k] = gammakx/gammaksum;
		//} else {
			//theta_arr[k] = 0.0;
		//}
	
		//TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		//TRYCXX( ISDestroy(&gammak_is) );
	//}	

	///* restore arrays */
	//TRYCXX( VecRestoreArray(theta_Vec,&theta_arr) );

	//TRYCXX( PetscBarrier(NULL));

	
	this->timer_solve.stop(); 

	LOG_FUNC_END
}

template<>
void EntropySolver<PetscVector>::compute_moments() {
	LOG_FUNC_BEGIN

	this->timer_compute_moments.start(); 

	Vec x_Vec = entropydata->get_x()->get_vector();
	Vec x_power_Vec = x_power->get_vector();
	Vec x_power_gammak_Vec = x_power_gammak->get_vector();
	Vec gamma_Vec = entropydata->get_gamma()->get_vector();
	Vec moments_Vec = moments->get_vector();
	
	Vec gammak_Vec;
	IS gammak_is;
	
	TRYCXX( VecCopy(x_Vec,x_power_Vec) ); /* x^1 */
	
	double *moments_arr, mysum, gammaksum;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );
	for(int km=0; km < entropydata->get_Km(); km++){
		
		for(int k=0;k<entropydata->get_K();k++){
			/* get gammak */
			this->entropydata->get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

			/* compute x_power_gammak */
			TRYCXX( VecPointwiseMult(x_power_gammak_Vec, gammak_Vec, x_power_Vec) ); /* x_power_gammak = x_power.*gammak */

			/* compute gammaksum */
			TRYCXX( VecSum(gammak_Vec, &gammaksum) );
			TRYCXX( VecSum(x_power_gammak_Vec, &mysum) );

			/* store computed moment */
			if(gammaksum != 0){
				moments_arr[k*this->entropydata->get_Km() + km] = mysum/gammaksum;
			} else {
				moments_arr[k*this->entropydata->get_Km() + km] = 0.0;
			}
	
			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );	
		}
		
		TRYCXX( VecPointwiseMult(x_power_Vec, x_Vec, x_power_Vec) ); /* x_power = x_power.*x */
	}
	
	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	this->timer_compute_moments.stop(); 
	
	LOG_FUNC_END
}


}
} /* end namespace */

#endif
