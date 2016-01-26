#ifndef PROJECTION_H
#define	PROJECTION_H

#include "common.h"

void get_projection(GammaVector<Scalar> *x, int K);
void get_projection(GammaVector<Scalar> *x, int K, double *time_to_add);

void get_projection_sub(Scalar *x_sub, int n);
void sort_bubble(Scalar *x, int n);

#endif

