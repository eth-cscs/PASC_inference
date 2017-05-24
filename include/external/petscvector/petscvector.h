#ifndef PASC_PETSCVECTOR_H
#define	PASC_PETSCVECTOR_H

#include "pascinference.h"

#include "petsc.h"
#include "petscsys.h"

#ifdef USE_CUDA
	 #include "petsccuda.h"
	 #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

/* the petscvector variant of general.h */
#include "external/petscvector/common/common.h"
#include "external/petscvector/algebra/algebra.h"
#include "external/petscvector/data/data.h"
#include "external/petscvector/model/model.h"
#include "external/petscvector/solver/solver.h"

#endif
