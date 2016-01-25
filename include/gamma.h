#ifndef GAMMA_H
#define	GAMMA_H

class Gamma;
class QPSolver;

#include "common.h"
#include "data.h"
#include "theta.h"
#include "qpsolver.h"

class Gamma {
		int K; /* number of gamma components*/
		int T; /* length of time serie */
				
	public:
		GammaVector<Scalar> gamma_vec; /* long vector with gamma [gamma_1, ... gamma_K] */ //TODO: should be private

		void init(Data, int); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print();
		void print(int nmb_of_spaces);

		void prepare_random();
		void prepare_uniform();		

		void compute(QPSolver *qpsolver, Data data, Theta theta);

		/* GET functions */
		int get_T();
		int get_K();

		// TODO: this should be somewhere else, Model?
		void compute_gk(GammaVector<Scalar> *g, Data *data, Theta *theta);
	
};

#endif

