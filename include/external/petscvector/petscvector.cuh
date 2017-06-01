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

#endif


