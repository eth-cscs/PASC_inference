#ifndef DATA_H
#define	DATA_H

#include "common.h"

class Data {
		int dim; /* number of data components */

		int global_size; /* global data size with overlap */
		int local_size; /* local data size with overlap */

		int local_begin, local_end; /* ownership range */

		PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	public:
		Vec *data_vecs; /* array with data vectors, TODO: should be private */

		void init(int dim, int global_size); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print(PetscViewer v);
		void get_covtrace(PetscScalar *covtrace);
		
		/* GET functions */
		int get_local_size();
		int get_global_size();
		int get_local_begin();
		int get_local_end();
		int get_dim();

	
};



#endif
