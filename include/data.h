#ifndef DATA_H
#define	DATA_H

#include "common.h"

class Data {
		PetscInt dim; /* number of data components */

		PetscInt global_size; /* global data size with overlap */
		PetscInt local_size; /* local data size with overlap */

		PetscInt local_begin, local_end; /* ownership range */

		PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	public:
		Vec *data_vecs; /* array with data vectors, TODO: should be private */

		PetscErrorCode init(PetscInt dim, PetscInt global_size); /* TODO: should be constructor */
		PetscErrorCode finalize(); /* TODO: should be destructor */

		PetscErrorCode print(PetscViewer v);
		PetscErrorCode get_covtrace(PetscScalar *covtrace);
		
		/* GET functions */
		PetscInt get_local_size();
		PetscInt get_global_size();
		PetscInt get_local_begin();
		PetscInt get_local_end();
		PetscInt get_dim();

	
};



#endif
