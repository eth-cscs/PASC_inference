#ifndef OPERATIONS_H
#define	OPERATIONS_H

#include "common.h"

Scalar get_dot(GammaVector<Scalar> x, GammaVector<Scalar> y);
Scalar get_dot(GammaVector<Scalar> x, GammaVector<Scalar> y, double *time_to_add);

void get_dot(Scalar *xx, GammaVector<Scalar> x, GammaVector<Scalar> y);
void get_dot(Scalar *xx, GammaVector<Scalar> x, GammaVector<Scalar> y, double *time_to_add);

void get_Ax_laplace(GammaVector<Scalar> *Ax, GammaVector<Scalar> x); 
void get_Ax_laplace(GammaVector<Scalar> *Ax, GammaVector<Scalar> x, double *time_to_add); 


#endif
