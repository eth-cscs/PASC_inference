#ifndef PASC_PETSCVECTOR_INITIALIZE_H
#define	PASC_PETSCVECTOR_INITIALIZE_H

#include "petsc.h"
#include "petscsys.h"
#ifdef USE_CUDA
 /* include additional petsc-cuda stuff */
 #include "petsccuda.h"
#endif

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/common/initialize.h"

namespace pascinference {
namespace common {

extern bool PETSC_INITIALIZED;

/* for loading PETSc options */
extern char **argv_petsc;
extern int argc_petsc;

template<> bool Initialize<PetscVector>(int argc, char *argv[]);
template<> void Finalize<PetscVector>();
template<> void allbarrier<PetscVector>();

}
} /* end of namespace */


#endif
