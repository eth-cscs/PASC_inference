#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"
#include "theta.h"
#include "data.h"
#include "operations.h"

class QPSolver {
	protected:
		Data *data;
		Gamma *gamma;
		Theta *theta;

		Scalar eps_sqr;
		int it;
		int hess_mult;
		
		/* data of QP problem */
		GammaVector<Scalar> *bs; /* rhs */
		GammaVector<Scalar> *gs; /* gradient */
		GammaVector<Scalar> *ds; /* projected gradient */
		GammaVector<Scalar> *Ads; /* A*ds */
		
		/* private functions */
		void project(GammaVector<Scalar> **x);
		void project_sub(GammaVector<Scalar> *x_sub);
		void sort_bubble(GammaVector<Scalar> *x);
		
		double time_projection; /* the sum of time necessary to perform projections */
		double time_matmult; /* the sum of time necessary to perform matrix multiplication */
		double time_dot; /* the sum of time necessary to compute dot_products */
		double time_update; /* total time of vector updates */
		double time_total; /* total time of SPG algorithm */
		
	public:
		QPSolver(Data*, Gamma *, Theta *, Scalar);
		void init();
		void finalize();

		void compute_b();
		void solve();
		Scalar get_function_value(GammaVector<Scalar> *x);
		Scalar get_function_value(GammaVector<Scalar> *x, bool use_gradient);

		void print();
		void print(int nmb_of_spaces);

		int get_T();
		int get_K();
		int get_dim();

		int get_it();
		int get_hessmult();
		double get_time_projection();
		double get_time_matmult();
		double get_time_dot();
		double get_time_update();
		double get_time_total();
		double get_time_other();

};



#endif
