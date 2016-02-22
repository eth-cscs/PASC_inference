#ifndef OPERATIONS_H
#define	OPERATIONS_H

#include "common.h"

namespace pascinference {

/* HostVector */
Scalar get_dot(minlin::threx::HostVector<Scalar> x, minlin::threx::HostVector<Scalar> y);
void get_dot(Scalar *xx, minlin::threx::HostVector<Scalar> x, minlin::threx::HostVector<Scalar> y);
void get_Ax_laplace(minlin::threx::HostVector<Scalar> &Ax, minlin::threx::HostVector<Scalar> x, int K); 
void get_Ax_laplace(minlin::threx::HostVector<Scalar> &Ax, minlin::threx::HostVector<Scalar> x, int K, Scalar alpha); 

/* DeviceVector */
#ifdef USE_GPU
	Scalar get_dot(minlin::threx::DeviceVector<Scalar> x, minlin::threx::DeviceVector<Scalar> y);
	void get_dot(Scalar *xx, minlin::threx::DeviceVector<Scalar> x, minlin::threx::DeviceVector<Scalar> y);
	void get_Ax_laplace(minlin::threx::DeviceVector<Scalar> &Ax, minlin::threx::DeviceVector<Scalar> x, int K); 
	void get_Ax_laplace(minlin::threx::DeviceVector<Scalar> &Ax, minlin::threx::DeviceVector<Scalar> x, int K, Scalar alpha); 
#endif

/* PetscVector */
#ifdef USE_PETSC
	Scalar get_dot(petscvector::PetscVector &x, petscvector::PetscVector &y);
	void get_dot(Scalar *xx, petscvector::PetscVector &x, petscvector::PetscVector &y);
	void get_Ax_laplace(petscvector::PetscVector &Ax, petscvector::PetscVector &x, int K); 
	void get_Ax_laplace(petscvector::PetscVector &Ax, petscvector::PetscVector &x, int K, Scalar alpha); 
#endif

} /* end of namespace */


#endif
