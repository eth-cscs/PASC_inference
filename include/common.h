/*! This class should include all common functions and includes */

#ifndef COMMON_H
#define	COMMON_H

#define datan 2
#define gammaK 3
#define max_s_steps 1000
#define deltaL_eps 0.00001

#define PRINT_DATA 0

/* define several solver type options */
typedef enum {
	QPSOLVER_PERMON,
	QPSOLVER_PROJECTIONSTEP
} QPSolverType;
static const char *const QPSolverType_names[] = {"PERMON","PROJECTIONSTEP", "QPSolverType","QPSOLVER_",0}; /* name of the options in console parameters */


/* include PETSc */
#include "petsc.h"
#include <petsctime.h> /* time management */

/* include common c++ header files */
#include <iostream>
//#include <complex>

using namespace std;

/* general utils */
void Initialize(int, char**);
void Finalize();


#endif
