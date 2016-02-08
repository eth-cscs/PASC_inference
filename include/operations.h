#ifndef OPERATIONS_H
#define	OPERATIONS_H

#include "common.h"

Scalar get_dot(GammaVector x, GammaVector y);
void get_dot(Scalar *xx, GammaVector x, GammaVector y);

void get_Ax_laplace(GammaVector& Ax, GammaVector x, int K); 
void get_Ax_laplace(GammaVector& Ax, GammaVector x, int K, Scalar alpha); 


#endif
