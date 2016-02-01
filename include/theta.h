#ifndef THETA_H
#define	THETA_H

class Theta;

#include "common.h"
#include "data.h"
#include "gamma.h"

class Theta {
		int dim; /* number of data components = n */
		int K; /* number of gamma components */

	public:
		ThetaVector<Scalar> theta_vec; /* vector with theta */

		void init(int dim, int K); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */
		void print();
		void print(int nmb_of_spaces);
		void compute(DataVector<Scalar> data_vec, Gamma gamma);
		
		
		/* GET functions */
		int get_dim();
		int get_K();

	
};

#endif
