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

		Gamma gamma; /* characteristic function of clusters */
		Theta theta; /* model parameters */

		Timer timer_gamma; /* for gamma manipulation */
		Timer timer_theta; /* for theta manipulation */

		
	public:

		void init(int dim, int T, int K); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */
		void print();
		void print_timers();

		void compute_theta(DataVector data_vec);
		void compute_gamma(DataVector data_vec);

		/* GET functions */
		double get_function_value();
		int get_dim();
		int get_T();
		int get_K();
		Gamma get_gamma(); // TODO: temp, just for test
		Theta get_theta(); // TODO: temp, just for test
	
};

class Model_kmeans : public Model {
	private:
	
	public:

	
};


#endif
