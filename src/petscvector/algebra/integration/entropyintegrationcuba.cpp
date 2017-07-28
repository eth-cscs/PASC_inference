#include "external/petscvector/algebra/integration/entropyintegrationcuba.h"

#ifdef USE_CUBA

namespace pascinference {
namespace algebra {

template<>
void EntropyIntegrationCuba<PetscVector>::compute(double *integrals_arr, double *lambda_arr, int Km_max) {
	LOG_FUNC_BEGIN

	timer.start();
	
	ExternalContent<PetscVector>::Integrator integrator(this->integration_type, this->xdim, this->number_of_integrals, this->integration_mineval, this->integration_maxeval, this->integration_nstart, this->integration_nincrease, this->debug_print_integration);

	/* setting to compute normalization constant */
	ExternalContent<PetscVector>::ExtraParameters xp(Mom, LM, D, 0.0);
	integrator.USERDATA = &xp;

	
		cubareal *computed_integrals;
	computed_integrals = integrator.compute();
	
	if(Km_max < 0){
		Km_max = this->number_of_moments;
	}
	
	column_vector lambda_Dlib(this->number_of_moments);

	/* from arr to Dlib-vec */
	for(int km=0;km<this->number_of_moments;km++){
		lambda_Dlib(km) = lambda_arr[km];
	}

	/* compute integrals */
    for(int km = 0; km<Km_max;km++){
		auto mom_function = [&](double x)->double { return this->externalcontent->gg(x, km, lambda_Dlib);};
		integrals_arr[km] = dlib::integrate_function_adapt_simp(mom_function, -1.0, 1.0, this->eps);
	}

	timer.stop();
	
	LOG_FUNC_END
}






/* ------------ ExtraParameters ----------------- */
EntropyIntegrationCuba<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters(){
    D = 0.0;
    eps= 0.0;
    Mom = 0.0;
    LM = 0.0;
}

EntropyIntegrationCuba<PetscVector>::ExternalContent::ExtraParameters::ExtraParameters(double *_Mom, double *_LM, double *_D, double _eps){
    Mom = _Mom;
    LM = _LM;
    eps = _eps;
    D = _D;
}

void EntropyIntegrationCuba<PetscVector>::ExternalContent::ExtraParameters::Copy(ExtraParameters &_ExtraParameters){
    Mom = _ExtraParameters.Mom;
    eps = _ExtraParameters.eps;
    D = _ExtraParameters.D;
    LM = _ExtraParameters.LM;
}

EntropyIntegrationCuba<PetscVector>::ExternalContent::ExtraParameters::~ExtraParameters(){
}


/* ------------ Integrator ----------------- */
int EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::Integrand(const int *ndim, const cubareal xx[],
                     const int *ncomp, cubareal ff2[], void *userdata) {

    ExtraParameters* xp = (ExtraParameters*)userdata;
    dlib::matrix<double> D = xp->D;
    long d = D.nc();
    long n = D.nr();

    column_vector LM = xp->LM;

    double V = 0.0;
    double p = 0.0;

    for (int i = 0; i < n; i++){
        p = 1.0;
        for (int j = 0; j < d; j++){
            p = p*pow(xx[j], D(i,j));
		}
        V = V - p*LM(i);
    }

	/* ff2[0] - type = 0 */
    ff2[0] = exp(V);

    /* ff2[1-n] - for gradient, type =2 */
	for(int order = 0; order < n; order++){
		p = 1.0;
		for (int j = 0; j < d; j++){
			p = p*pow(xx[j], D(order,j));
		}
        ff2[1+order] = p*exp(V);
    }

	/* ff2[n+1 - n+1+n*(n+1)/2] - for Hessian, type = 3 */
	int counter = 1+n;
	for(int order = 0; order < n; order++){
		for(int order2 = order; order2 < n; order2++){
			p = 1.0;
			row_vector t = rowm(D,order) + rowm(D,order2);
			for (int j = 0; j < d; j++){
				p = p*pow(xx[j], t(j));
			}
			ff2[counter] = p*exp(V);
			counter++;
		}
	}

    return 0;
}

EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::Integrator(int integration_type, int ndim, int ncomp, int integration_mineval, int integration_maxeval, int integration_nstart, int integration_nincrease, bool debug_print_integration){
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

	this->integration_type = integration_type;
	this->integral = new cubareal[ncomp];
	this->error = new cubareal[ncomp];
	this->prob = new cubareal[ncomp];

	this->debug_print_integration = debug_print_integration;

}

cubareal* EntropySolverDlib<PetscVector>::ExternalContent::Integrator::compute() {
	switch(this->integration_type){
		case 0: this->computeVegas(); break;
		case 1: this->computeSuave(); break;
		case 2: this->computeDivonne(); break;
		case 3: this->computeCuhre(); break;
	}
    return integral;
}

void EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::computeVegas(){
	timer.start();

    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
          EPSREL, EPSABS, VERBOSE, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, SPIN,
          &neval, &fail, integral, error, prob);

	timer.stop();

    if(this->debug_print_integration){
		coutMaster << "VEGAS RESULT: " << neval << " neval," << fail << " fail" << std::endl;
		for(int comp = 0; comp < NCOMP; comp++ ){
			coutMaster.push();
            coutMaster << "value: " << (double)integral[comp] << ", error: " << (double)error[comp] << ", prob: " << (double)prob[comp] << std::endl;
			coutMaster.pop();
		}
	}
}

void EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::computeSuave(){
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

void EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::computeDivonne(){
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

void EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::computeCuhre(){
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

EntropyIntegrationCuba<PetscVector>::ExternalContent::Integrator::~Integrator(){
	free(this->integral);
	free(this->error);
	free(this->prob);
}



}
}

#endif
