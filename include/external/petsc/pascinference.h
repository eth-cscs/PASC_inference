#ifndef PASCINFERENCE_PETSC_H
#define	PASCINFERENCE_PETSC_H

#include "pascinference.h"

#include "petsc.h"
#include "petscsys.h"

#ifdef USE_CUDA
/* anselm hotfix */
	typedef struct _p_PetscCUDAIndices* PetscCUDAIndices;
	typedef struct _p_VecScatterCUDAIndices_StoS* VecScatterCUDAIndices_StoS;
	typedef struct _p_VecScatterCUDAIndices_PtoP* VecScatterCUDAIndices_PtoP;
	PETSC_EXTERN PetscErrorCode VecCUDACopyToGPUSome_Public(Vec,PetscCUDAIndices);
	PETSC_EXTERN PetscErrorCode VecCUDACopyFromGPUSome_Public(Vec,PetscCUDAIndices);
	PETSC_EXTERN PetscErrorCode VecScatterInitializeForGPU(VecScatter,Vec,ScatterMode);
	PETSC_EXTERN PetscErrorCode VecScatterFinalizeForGPU(VecScatter);
	PETSC_EXTERN PetscErrorCode VecCreateSeqCUDA(MPI_Comm,PetscInt,Vec*);
	PETSC_EXTERN PetscErrorCode VecCreateMPICUDA(MPI_Comm,PetscInt,PetscInt,Vec*);
#endif

/* PetscVector */
#include "external/petsc/algebra/vector/petscvector.h"

#endif
