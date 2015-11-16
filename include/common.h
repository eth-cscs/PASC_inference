/*! This class should include all common functions and includes */

#ifndef COMMON_H
#define	COMMON_H

#define dataN 1000
#define datan 2
#define gammaK 3
#define max_s_steps 1000
#define deltaL_eps 0.0001

#define PRINT_DATA 0

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
