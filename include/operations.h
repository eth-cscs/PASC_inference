#ifndef OPERATIONS_H
#define	OPERATIONS_H

#include "common.h"

/* HostVector */
Scalar get_dot(HostVector<Scalar> x, HostVector<Scalar> y);
void get_dot(Scalar *xx, HostVector<Scalar> x, HostVector<Scalar> y);
void get_Ax_laplace(HostVector<Scalar> &Ax, HostVector<Scalar> x, int K); 
void get_Ax_laplace(HostVector<Scalar> &Ax, HostVector<Scalar> x, int K, Scalar alpha); 

/* DeviceVector */
#ifdef USE_GPU
	Scalar get_dot(DeviceVector<Scalar> x, DeviceVector<Scalar> y);
	void get_dot(Scalar *xx, DeviceVector<Scalar> x, DeviceVector<Scalar> y);
	void get_Ax_laplace(DeviceVector<Scalar> &Ax, DeviceVector<Scalar> x, int K); 
	void get_Ax_laplace(DeviceVector<Scalar> &Ax, DeviceVector<Scalar> x, int K, Scalar alpha); 
#endif

/* PetscVector */
#ifdef USE_PETSC
	Scalar get_dot(PetscVector &x, PetscVector &y);
	void get_dot(Scalar *xx, PetscVector &x, PetscVector &y);
	void get_Ax_laplace(PetscVector &Ax, PetscVector &x, int K); 
	void get_Ax_laplace(PetscVector &Ax, PetscVector &x, int K, Scalar alpha); 
#endif


#endif
