/*! This class should include all common functions and includes */

#ifndef COMMON_H
#define	COMMON_H

#define dataN 10
#define datan 2
#define gammaK 3
#define eps_squared 10

#define PRINT_DATA 1

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
