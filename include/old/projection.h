#ifndef PROJECTION_H
#define	PROJECTION_H

#include "common.h"

namespace pascinference {
	
/* -------- HostVector ---------- */
void get_projection(minlin::threx::HostVector<Scalar> & x, int K);

void get_projection_sub(Scalar *x_sub, int n);
void sort_bubble(Scalar *x, int n);


/* -------- DeviceVector ---------- */
#ifdef USE_GPU
	void get_projection(minlin::threx::DeviceVector<Scalar> & x, int K);

	__device__
	void device_sort_bubble(Scalar *x, int n);

	__global__
	void kernel_get_projection_sub(Scalar *x, int T, const int K);
#endif

/* -------- PetscVector ---------- */
#ifdef USE_PETSC
	void get_projection(petscvector::PetscVector & x, int K);
#endif

} /* end of namespace */



#endif

