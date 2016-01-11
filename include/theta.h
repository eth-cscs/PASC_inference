#ifndef THETA_H
#define	THETA_H

class Theta;

#include "common.h"
#include "data.h"
#include "gamma.h"

class Theta {
		int dim_data; /* number of data components = n */
		int dim_gamma; /* number of gamma components = K */

	public:
		float *theta_arr; /* array with theta, should be same in every proc */

		void init(Data data, Gamma gamma); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */
		void print(PetscViewer v);
		void compute(Data data, Gamma gamma);
		
		
		/* GET functions */
		int get_dim_data();
		int get_dim_gamma();

	
};

#endif
