#ifndef GAMMA_H
#define	GAMMA_H

class Gamma;
class QPSolver;

#include "common.h"
#include "data.h"
#include "theta.h"
#include "qpsolver.h"

class Gamma {
		int dim;
		int K; /* number of gamma components*/
		int T; /* length of time serie */

		QPSolver qpsolver;
				
	public:
		GammaVector gamma_vec; /* long vector with gamma [gamma_1, ... gamma_K] */ //TODO: should be private

		void init(int dim, int T, int K); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print();
		void print_timers();
		void print(int nmb_of_spaces);

		void prepare_random();
		void prepare_uniform();		

		void compute(DataVector data_vec, Theta theta);

		/* GET functions */
		int get_T();
		int get_K();
		int get_function_value();
		GammaVector get_gamma_vec();

		// TODO: this should be somewhere else, Model?
		void compute_gk(GammaVector& g, DataVector data_vec, Theta theta);
	
};

#endif

