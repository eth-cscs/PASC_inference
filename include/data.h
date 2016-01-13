#ifndef DATA_H
#define	DATA_H

#include "common.h"

class Data {
		int dim; /* number of data components */
		int T; /* length of time serie */

	public:
		DataVector<Scalar> *data_vecs; /* array with data vectors, TODO: should be private */

		void init(int dim, int T); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print();
//		void get_covtrace(PetscScalar *covtrace);
		
		/* GET functions */
		int get_dim();
		int get_T();
	
};



#endif
