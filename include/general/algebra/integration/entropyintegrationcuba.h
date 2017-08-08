/** @file entropyintegrationcuba.h
 *  @brief Computes numerical integrals using cuba library
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_ENTROPYINTEGRATIONCUBA_H
#define	PASC_ENTROPYINTEGRATIONCUBA_H

#ifdef USE_CUBA
/* if we are not using CUBA, then this class does not make any sence */

#include <string>
#include <iostream>
#include "general/algebra/integration/entropyintegration.h"
#include "general/common/consoleoutput.h"
#include "general/common/logging.h"
#include "general/common/timer.h"

#define ENTROPYINTEGRATIONCUBA_DEFAULT_TYPE 0
#define ENTROPYINTEGRATIONCUBA_DEFAULT_MINEVAL 0
#define ENTROPYINTEGRATIONCUBA_DEFAULT_MAXEVAL 5e4
#define ENTROPYINTEGRATIONCUBA_DEFAULT_NSTART 1e4
#define ENTROPYINTEGRATIONCUBA_DEFAULT_NINCREASE 1e4

#define ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION false

/* define number of accellerators used for integration */ 
#ifdef USE_GPU
	#define ENTROPYINTEGRATIONCUBA_DEFAULT_CUDAACCEL 0
#else
	#define ENTROPYINTEGRATIONCUBA_DEFAULT_CUDAACCEL 1
#endif
#define ENTROPYINTEGRATIONCUBA_DEFAULT_CUDAACCELMAX 1000

/* include Cuba stuff */
#include "cuba.h"

namespace pascinference {
using namespace common;

namespace algebra {

template<class VectorBase>
class EntropyIntegrationCuba : public EntropyIntegration<VectorBase> {
	public:
		class ExternalContent;

		class Integrator {
			public:
				int NDIM; 		/**< number of dimensions of integral */
				int NCOMP;		/**< number of components of the integrand */
				int NVEC;
				double EPSREL;	/**< requested relative accuracy */
				double EPSABS;	/**< requested absolute accuracy */
				int VERBOSE; //log output
				int LAST;
				int SEED;
				int MINEVAL;	/**< minimum number of integrand evaluations */
				int MAXEVAL;	/**< maximum number of integrand evaluations allowed */
				int NSTART;
				int NINCREASE;
				int NBATCH;
				int GRIDNO;
				char* STATEFILE = NULL;
				void* SPIN = NULL;
				int NNEW;
				int NMIN;
				double FLATNESS;
				void* USERDATA = NULL; //this is to pass extra parameters to integral

				int KEY1;
				int KEY2;
				int KEY3;
				int MAXPASS;
				double BORDER;
				double MAXCHISQ;
				double MINDEVIATION;
				int NGIVEN;
				int LDXGIVEN;
				int NEXTRA;
				int KEY;

				int CUBAACCEL;
				int CUBAACCELMAX;

				int integration_type;

				int comp, nregions, neval, fail;

				cubareal *integral;
				cubareal *error;
				cubareal *prob;

				bool debug_print_integration;

				Integrator(int integration_type, int ndim, int ncomp, int integration_mineval, int integration_maxeval, int integration_nstart, int integration_nincrease, int integration_cubaaccel, int integration_cubaaccelmax, bool debug_print_integration = ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION);
				~Integrator();

				//four methods of integration implemented in CUBA library,
				//more info at http://www.feynarts.de/cuba/
				void computeVegas();
				void computeSuave();
				void computeDivonne();
				void computeCuhre();
				cubareal *compute();

				static int Integrand(const int *ndim, const double xx[], const int *ncomp, cubareal ff2[], void *userdata);
		};

		class ExtraParameters {
			public:
				double *LM;
				double eps;
				int *matrix_D_arr;
				int xdim; /* number of columns in D */
				int number_of_moments; /* number of rows in D */

				ExtraParameters();
				ExtraParameters(double *_LM, int *_matrix_D_arr, int _xdim, int _number_of_moments, double _eps);
				~ExtraParameters();
				void Copy(ExtraParameters& _ExtraParameters);
		};

		int type;					/**< integration type [0=Vegas,1=Suave,2=Divonne,3=Cuhre] */
		int mineval;				/**< the minimum number of integrand evaluations */
		int maxeval;				/**< the maximum number of integrand evaluations */
		int nstart;					/**< number of integrand evaluations to start with */
		int nincrease;				/**< the increase in number of integrand evaluations */

		int cudaaccel;
		int cudaaccelmax;

		bool debug_print_integration;

		void set_settings_from_console();


	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

	public:
		EntropyIntegrationCuba(EntropyData<VectorBase> *entropydata, double new_eps);
		~EntropyIntegrationCuba();

		virtual std::string get_name() const;
		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void compute(double *integrals_out, double *lambda, int Km_max = -1);

		std::string get_integration_type_name(int integration_type) const;

};

}
} /* end of namespace */


/* ----- IMPLEMENTATION ----- */
namespace pascinference {
namespace algebra {

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::set_settings_from_console() {

	consoleArg.set_option_value("entropyintegrationcuba_type", &this->type, ENTROPYINTEGRATIONCUBA_DEFAULT_TYPE);
	consoleArg.set_option_value("entropyintegrationcuba_mineval", &this->mineval, ENTROPYINTEGRATIONCUBA_DEFAULT_MINEVAL);
	consoleArg.set_option_value("entropyintegrationcuba_maxeval", &this->maxeval, ENTROPYINTEGRATIONCUBA_DEFAULT_MAXEVAL);
	consoleArg.set_option_value("entropyintegrationcuba_nstart", &this->nstart, ENTROPYINTEGRATIONCUBA_DEFAULT_NSTART);
	consoleArg.set_option_value("entropyintegrationcuba_nincrease", &this->nincrease, ENTROPYINTEGRATIONCUBA_DEFAULT_NINCREASE);

	consoleArg.set_option_value("entropyintegrationcuba_cudaaccel", &this->cudaaccel, ENTROPYINTEGRATIONCUBA_DEFAULT_CUDAACCEL);
	consoleArg.set_option_value("entropyintegrationcuba_cudaaccelmax", &this->cudaaccelmax, ENTROPYINTEGRATIONCUBA_DEFAULT_CUDAACCELMAX);

	consoleArg.set_option_value("entropyintegrationcuba_debug_print_integration",&debug_print_integration, ENTROPYINTEGRATIONCUBA_DEFAULT_DEBUG_PRINT_INTEGRATION);
}

/* constructor */
template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::EntropyIntegrationCuba(EntropyData<VectorBase> *entropydata, double new_eps) : EntropyIntegration<VectorBase>(entropydata, new_eps) {
	LOG_FUNC_BEGIN

	/* load parameters from console */
	set_settings_from_console();

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::~EntropyIntegrationCuba(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EntropyIntegrationCuba<VectorBase>::get_name() const {
	std::string return_value = "EntropyIntegrationCuba<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
std::string EntropyIntegrationCuba<VectorBase>::get_integration_type_name(int integration_type) const {
	std::string return_string = "undefined";
	switch(integration_type){
		case 0: return_string = "Vegas"; break;
		case 1: return_string = "Suave"; break;
		case 2: return_string = "Divonne"; break;
		case 3: return_string = "Cuhre"; break;
	}
	return return_string;
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;

	output <<  " - number of moments        : " << this->entropydata->get_number_of_moments() << std::endl;
	output <<  " - xdim                     : " << this->entropydata->get_xdim() << std::endl;
	output <<  " - eps                      : " << this->eps << std::endl;

	output <<  " - integration_type         : " << get_integration_type_name(this->type) << std::endl;
	output <<  " - integration_mineval      : " << this->mineval << std::endl;
	output <<  " - integration_maxeval      : " << this->maxeval << std::endl;
	output <<  " - integration_nstart       : " << this->nstart << std::endl;
	output <<  " - integration_nincrease    : " << this->nincrease << std::endl;

	output <<  " - integration_cubaaccel    : " << this->cudaaccel << std::endl;
	output <<  " - integration_cubaaccelmax : " << this->cudaaccelmax << std::endl;

	output.synchronize();

	LOG_FUNC_END
}

/* print info about integration */
template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	output_global <<  " - number of moments       : " << this->entropydata->get_number_of_moments() << std::endl;
	output_global <<  " - xdim                    : " << this->entropydata->get_xdim() << std::endl;
	output_global <<  " - eps                     : " << this->eps << std::endl;

	output_global <<  " - integration_type        : " << get_integration_type_name(this->type) << std::endl;
	output_global <<  " - integration_mineval     : " << this->mineval << std::endl;
	output_global <<  " - integration_maxeval     : " << this->maxeval << std::endl;
	output_global <<  " - integration_nstart      : " << this->nstart << std::endl;
	output_global <<  " - integration_nincrease   : " << this->nincrease << std::endl;

	output_global <<  " - integration_cubaaccel    : " << this->cudaaccel << std::endl;
	output_global <<  " - integration_cubaaccelmax : " << this->cudaaccelmax << std::endl;

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::compute(double *integrals_out, double *lambda, int Km_max) {
	LOG_FUNC_BEGIN

	this->timer.start();

	Integrator integrator(this->type, this->entropydata->get_xdim(), Km_max, this->mineval, this->maxeval, this->nstart, this->nincrease, this->cudaaccel, this->cudaaccelmax, this->debug_print_integration);

	/* setting to compute normalization constant */
	ExtraParameters xp(lambda, this->entropydata->get_matrix_D(), this->entropydata->get_xdim(), this->entropydata->get_number_of_moments()-1, 0.0);
	integrator.USERDATA = &xp;

	cubareal *computed_integrals;
	computed_integrals = integrator.compute();

	for(int i=0;i<Km_max;i++){
		integrals_out[i] = computed_integrals[i];
	}

    if(debug_print_integration){
        coutAll << "lambda:    " << print_array(lambda, this->entropydata->get_number_of_moments()-1) << std::endl;
        coutAll << "integrals: " << print_array(computed_integrals, Km_max) << std::endl;
        coutAll.synchronize();
    }

	this->timer.stop();

	LOG_FUNC_END
}





/* ------------ Integrator ----------------- */
template<class VectorBase>
int EntropyIntegrationCuba<VectorBase>::Integrator::Integrand(const int *ndim, const double xx[],
                     const int *ncomp, cubareal ff2[], void *userdata) {

    ExtraParameters* xp = (ExtraParameters*)userdata;
    int *D = xp->matrix_D_arr;
    long d = xp->xdim;
    long n = xp->number_of_moments;

    double *LM = xp->LM;

    double V = 0.0;
    double p = 0.0;

    for (int i = 0; i < n; i++){
        p = 1.0;
        for (int j = 0; j < d; j++){
            p = p*pow(xx[j], D[(i+1)*d+j]);
		}

        V = V - p*LM[i];
    }

	/* ff2[0] - type = 0 */
    ff2[0] = exp(V);

    /* ff2[1-n] - for gradient, type =2 */
	for(int order = 0; order < n; order++){
		p = 1.0;
		for (int j = 0; j < d; j++){
			p = p*pow(xx[j], D[(order+1)*d+j]);
		}
        ff2[1+order] = p*ff2[0];
    }

	/* ff2[n+1 - n+1+n*(n+1)/2] - for Hessian, type = 3 */
	int counter = 1+n;
	for(int order = 0; order < n; order++){
		for(int order2 = order; order2 < n; order2++){
			p = 1.0;
			for(int j=0; j<d;j++){
				p = p*pow(xx[j], D[(order+1)*d+j] + D[(order2+1)*d+j]);
			}
			ff2[counter] = p*ff2[0];
			counter++;
		}
	}

    return 0;
}

template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::Integrator::Integrator(int integration_type, int ndim, int ncomp, int integration_mineval, int integration_maxeval, int integration_nstart, int integration_nincrease, int integration_cubaaccel, int integration_cubaaccelmax, bool debug_print_integration){
	//all this paramterers are from example file demo-c.c
	NDIM = ndim;
	NCOMP = ncomp;
	NVEC = 1;
	EPSREL = 1e-12;//1e-3;
	EPSABS = 1e-12;
	VERBOSE = 0; //log output
	LAST = 4;
	SEED = 0;
	MINEVAL = integration_mineval;
	MAXEVAL = integration_maxeval;//50000;
	NSTART = integration_nstart;//1000;
	NINCREASE = integration_nincrease;//500;
	NBATCH = 1000;
	GRIDNO = 0;
	STATEFILE = NULL;
	SPIN = NULL;
	NNEW = 1000;
	NMIN = 2;
	FLATNESS = 25.0;
	USERDATA = NULL; //this is to pass extra parameters to integral

	KEY1 = 47;
	KEY2 = 1;
	KEY3 = 1;
	MAXPASS = 5;
	BORDER = 0.0;
	MAXCHISQ = 10.0;
	MINDEVIATION = 0.25;
	NGIVEN = 0;
	LDXGIVEN = NDIM;
	NEXTRA = 0;

	KEY = 0;

	CUBAACCEL = integration_cubaaccel;
	CUBAACCELMAX = integration_cubaaccelmax;

	this->integration_type = integration_type;
	this->integral = new cubareal[ncomp];
	this->error = new cubareal[ncomp];
	this->prob = new cubareal[ncomp];

	this->debug_print_integration = debug_print_integration;

}

template<class VectorBase>
cubareal* EntropyIntegrationCuba<VectorBase>::Integrator::compute() {
	switch(this->integration_type){
		case 0: this->computeVegas(); break;
		case 1: this->computeSuave(); break;
		case 2: this->computeDivonne(); break;
		case 3: this->computeCuhre(); break;
	}
    return integral;
}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::Integrator::computeVegas(){
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
          EPSREL, EPSABS, VERBOSE, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, SPIN,
          &neval, &fail, integral, error, prob);

    if(this->debug_print_integration){
		coutMaster << "VEGAS RESULT: " << neval << " neval," << fail << " fail" << std::endl;
		for(int comp = 0; comp < NCOMP; comp++ ){
			coutMaster.push();
            coutMaster << "value: " << (double)integral[comp] << ", error: " << (double)error[comp] << ", prob: " << (double)prob[comp] << std::endl;
			coutMaster.pop();
		}
	}
}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::Integrator::computeSuave(){
    Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST, SEED,
              MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);

    if(this->debug_print_integration){
		coutMaster << "SUAVE RESULT: " << nregions << " nregions, " << neval << " neval," << fail << " fail" << std::endl;
		for(int comp = 0; comp < NCOMP; comp++ ){
			coutMaster.push();
            coutMaster << "value: " << (double)integral[comp] << ", error: " << (double)error[comp] << ", prob: " << (double)prob[comp] << std::endl;
			coutMaster.pop();
		}
	}

}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::Integrator::computeDivonne(){
    Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                BORDER, MAXCHISQ, MINDEVIATION,
                NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

    if(this->debug_print_integration){
		coutMaster << "DIVONNE RESULT: " << nregions << " nregions, " << neval << " neval," << fail << " fail" << std::endl;
		for(int comp = 0; comp < NCOMP; comp++ ){
			coutMaster.push();
            coutMaster << "value: " << (double)integral[comp] << ", error: " << (double)error[comp] << ", prob: " << (double)prob[comp] << std::endl;
			coutMaster.pop();
		}
	}

}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::Integrator::computeCuhre(){
    Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST,
              MINEVAL, MAXEVAL, KEY,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);

    if(this->debug_print_integration){
		coutMaster << "CUHRE RESULT: " << nregions << " nregions, " << neval << " neval," << fail << " fail" << std::endl;
		for(int comp = 0; comp < NCOMP; comp++ ){
			coutMaster.push();
            coutMaster << "value: " << (double)integral[comp] << ", error: " << (double)error[comp] << ", prob: " << (double)prob[comp] << std::endl;
			coutMaster.pop();
		}
	}

}

template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::Integrator::~Integrator(){
	free(this->integral);
	free(this->error);
	free(this->prob);
}


/* ------------ ExtraParameters ----------------- */
template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::ExtraParameters::ExtraParameters(){
    this->D = NULL;
    this->xdim = 1;
    this->number_of_moments = 0;
    eps= 0.0;
    LM = NULL;
}

template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::ExtraParameters::ExtraParameters(double *_LM, int *_matrix_D_arr, int _xdim, int _number_of_moments, double _eps){
    this->LM = _LM;
    this->eps = _eps;
    this->matrix_D_arr = _matrix_D_arr;
    this->xdim = _xdim;
    this->number_of_moments = _number_of_moments;
}

template<class VectorBase>
void EntropyIntegrationCuba<VectorBase>::ExtraParameters::Copy(ExtraParameters &_ExtraParameters){
    eps = _ExtraParameters.eps;
    matrix_D_arr = _ExtraParameters.matrix_D_arr;
    LM = _ExtraParameters.LM;
    xdim = _ExtraParameters.xdim;
    number_of_moments = _ExtraParameters.number_of_moments;
}

template<class VectorBase>
EntropyIntegrationCuba<VectorBase>::ExtraParameters::~ExtraParameters(){
}


}
} /* end namespace */


#endif /* USE_CUBA */

#endif



