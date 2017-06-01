#ifndef PASC_PETSCVECTOR_CUDA_H
#define	PASC_PETSCVECTOR_CUDA_H

#include <stdio.h> /* printf in cuda */

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <device_functions.h>

/* cuda error check */ 
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"\n\x1B[31mCUDA error:\x1B[0m %s %s \x1B[33m%d\x1B[0m\n\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

/* anselm hotfix */
/*
typedef struct _p_PetscCUDAIndices* PetscCUDAIndices;
typedef struct _p_VecScatterCUDAIndices_StoS* VecScatterCUDAIndices_StoS;
typedef struct _p_VecScatterCUDAIndices_PtoP* VecScatterCUDAIndices_PtoP;
PETSC_EXTERN PetscErrorCode VecCUDACopyToGPUSome_Public(Vec,PetscCUDAIndices);
PETSC_EXTERN PetscErrorCode VecCUDACopyFromGPUSome_Public(Vec,PetscCUDAIndices);
PETSC_EXTERN PetscErrorCode VecScatterInitializeForGPU(VecScatter,Vec,ScatterMode);
PETSC_EXTERN PetscErrorCode VecScatterFinalizeForGPU(VecScatter);
PETSC_EXTERN PetscErrorCode VecCreateSeqCUDA(MPI_Comm,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPICUDA(MPI_Comm,PetscInt,PetscInt,Vec*);
*/

#endif


