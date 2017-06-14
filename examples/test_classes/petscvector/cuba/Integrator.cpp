//
//  Integrator.cpp
//  Theta_problem_MV
//
//  Created by Ganna Marchenko on 08/06/17.
//
//

#include "Integrator.h"
#include "cuba.h"
#include "ExtraParameters.h"
#include "dlib/matrix.h"
using namespace dlib;

int Integrator::Integrand(const int *ndim, const cubareal xx[],
                     const int *ncomp, cubareal ff2[], void *userdata) {
    

    //double p1 = xx[0];
    //double p2 = xx[1];
    //double p3 = xx[2];
    ExtraParameters* xp = (ExtraParameters*)userdata;
    matrix<double> D = xp->D;
    long d = D.nc();
    long n = D.nr();
    int type = xp->type;
    column_vector LM = xp->LM;
    
    double V = 0.0;
    double p = 0.0;
    
    for (int i = 0; i < n; i++)
    {
        p = 1.0;
        for (int j = 0; j < d; j++)
            p = p*pow(xx[j], D(i,j));
        V = V - p*LM(i);
    }
    
    if (type == 0) //just exp density
        ff2[0] = exp(V);
    else if (type == 1)
        ff2[0] = exp(-xp->L0 + V);
    else if (type == 2) /* for gradient */
    {
        p = 1.0;
        for (int j = 0; j < d; j++)
            p = p*pow(xx[j], D(xp->order,j));
        ff2[0] = p*exp(V);
    }
    else if (type == 3) /* for Hessian */
    {
        p = 1.0;
        row_vector t = rowm(D,xp->order) + rowm(D,xp->order2);
        for (int j = 0; j < d; j++)
            p = p*pow(xx[j], t(j));
        ff2[0] = p*exp(V);
    }
    return 0;
}

Integrator::Integrator()
{
    //all this paramterers are from example file demo-c.c
    NDIM = 2; //dimensions of integral
    NCOMP = 1;
    NVEC = 1;
    EPSREL = 1e-3;
    EPSABS = 1e-12;
    VERBOSE = 0; //log output
    LAST = 4;
    SEED = 0;
    MINEVAL = 0;
    MAXEVAL = 50000;
    NSTART = 1000;
    NINCREASE = 500;
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

}

double Integrator::computeVegas()
{
    //printf("-------------------- Vegas test --------------------\n");
    
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
          EPSREL, EPSABS, VERBOSE, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, SPIN,
          &neval, &fail, integral, error, prob);
    
    //printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    //       neval, fail);
    //for( comp = 0; comp < NCOMP; ++comp )
    //    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //           (double)integral[comp], (double)error[comp], (double)prob[comp]);
    return (double)integral[0];
}

double Integrator::computeSuave()
{
    printf("\n-------------------- Suave test --------------------\n");
    
    Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST, SEED,
              MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);
    
    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
               nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
               (double)integral[comp], (double)error[comp], (double)prob[comp]);
    return (double)integral[0];
}

double Integrator::computeDivonne()
{
    printf("\n------------------- Divonne test -------------------\n");
    
    Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                BORDER, MAXCHISQ, MINDEVIATION,
                NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);
    
    printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
            nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
              (double)integral[comp], (double)error[comp], (double)prob[comp]);
    return (double)integral[0];
}

double Integrator::computeCuhre()
{
    printf("\n-------------------- Cuhre test --------------------\n");
    
    Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
              EPSREL, EPSABS, VERBOSE | LAST,
              MINEVAL, MAXEVAL, KEY,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);
    
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
               nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
               (double)integral[comp], (double)error[comp], (double)prob[comp]);
    return (double)integral[0];
}

Integrator::~Integrator()
{
}
