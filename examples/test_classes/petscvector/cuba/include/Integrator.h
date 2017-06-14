//
//  Integrator.h
//  Theta_problem_MV
//
//  Created by Ganna Marchenko on 08/06/17.
//
//

#ifndef __Theta_problem_MV__Integrator__
#define __Theta_problem_MV__Integrator__

#include <stdio.h>
#include "cuba.h"

class Integrator
{
public:
//all this paramterers are from example file demo-c.c
    int NDIM; //dimensions of integral
    int NCOMP;
    int NVEC;
    double EPSREL;
    double EPSABS;
    int VERBOSE; //log output
    int LAST;
    int SEED;
    int MINEVAL;
    int MAXEVAL;
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

    
    int comp, nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    
    Integrator();
    ~Integrator();
    
    //four methods of integration implemented in CUBA library,
    //more info at http://www.feynarts.de/cuba/
    double computeVegas();
    double computeSuave();
    double computeDivonne();
    double computeCuhre();
    
    static int Integrand(const int *ndim, const cubareal xx[],
                         const int *ncomp, cubareal ff2[], void *userdata);
    
};
#endif /* defined(__Theta_problem_MV__Integrator__) */
