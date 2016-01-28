#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"

#include "operations.h"
#include "projection.h"

class QPSolver {
	protected:
		Data *data;
		Gamma *gamma;
		Theta *theta;

		Scalar eps_sqr;

		int it;
		int it_all;
		int hessmult;
		int hessmult_all;
		
		/* data of QP problem */
		GammaVector<Scalar> g; /* gradient */
		GammaVector<Scalar> d; /* projected gradient */
		GammaVector<Scalar> Ad; /* A*ds */
		
		double time_projection; /* the sum of time necessary to perform projections */
		double time_matmult; /* the sum of time necessary to perform matrix multiplication */
		double time_dot; /* the sum of time necessary to compute dot_products */
		double time_update; /* total time of vector updates */
		double time_stepsize; /* total time of step-size computation */
		double time_fs; /* total time of manipulation with fs vector during iterations */
		double time_init; /* total time of initializing qpsolver */
		double time_total; /* total time of SPG algorithm */
		
	public:
		GammaVector<Scalar> b; /* rhs */ // TODO: this should be private

		QPSolver(Gamma *, Scalar);
		void init();
		void finalize();

		void solve();
		Scalar get_function_value();
		Scalar get_function_value(GammaVector<Scalar> x);
		Scalar get_function_value(GammaVector<Scalar> x, bool use_gradient);

		void print();
		void print(int nmb_of_spaces);

		int get_T();
		int get_K();
		int get_dim();

		int get_it();
		int get_it_all();
		int get_hessmult();
		int get_hessmult_all();
		double get_time_projection();
		double get_time_matmult();
		double get_time_dot();
		double get_time_update();
		double get_time_init();
		double get_time_stepsize();
		double get_time_fs();
		double get_time_total();
		double get_time_other();

};



#endif
