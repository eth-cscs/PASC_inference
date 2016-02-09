#ifndef PROJECTION_H
#define	PROJECTION_H

#include "common.h"

/* -------- DeviceVector ---------- */

void get_projection(DeviceVector<Scalar> & x, int K);

__device__
void device_sort_bubble(Scalar *x, int n);

__global__
void kernel_get_projection_sub(Scalar *x, int T, const int K);

/* -------- HostVector ---------- */
void get_projection(HostVector<Scalar> & x, int K);

void get_projection_sub(Scalar *x_sub, int n);
void sort_bubble(Scalar *x, int n);

/* -------- PetscVector ---------- */
void get_projection(PetscVector & x, int K);


#endif

