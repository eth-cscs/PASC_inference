#ifndef PROJECTION_H
#define	PROJECTION_H

#include "common.h"

void get_projection(GammaVector<Scalar> **x, int K);
void get_projection(GammaVector<Scalar> **x, int K, double *time_to_add);

void get_projection_sub(GammaVector<Scalar> *x_sub);
void sort_bubble(GammaVector<Scalar> *x);


#endif

