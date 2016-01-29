#ifndef MODEL_H
#define	MODEL_H

#include "common.h"
#include "gamma.h"
#include "theta.h"

class Model {
	protected:
		int dim; /* number of data components */
		int T; /* length of time-series */
		int K; /* number of clusters */

		Theta theta; /* model parameters */
		Gamma gamma; /* characteristic function of clusters */
		
	public:

		void init(int dim, int T, int K); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		/* GET functions */
		Gamma get_gamma(); // TODO: temp, just for test
		Theta get_theta(); // TODO: temp, just for test
	
};

class Model_kmeans : public Model {
	private:
	
	public:

	
};


#endif
