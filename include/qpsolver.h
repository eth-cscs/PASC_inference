#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"
#include "theta.h"
#include "data.h"

/* thrust tools */
#include <thrust/sort.h>

class QPSolver {
	protected:
		Data *data;
		Gamma *gamma;
		Theta *theta;

		Scalar eps_sqr;
		
		/* data of QP problem */
		GammaMatrix<Scalar> A_sub; /* Hessian */
		GammaVector<Scalar> *bs; /* rhs */
		GammaVector<Scalar> *gs; /* gradient */
		GammaVector<Scalar> *ds; /* projected gradient */
		GammaVector<Scalar> *Ads; /* A*ds */
		
		/* private functions */
		void project(GammaVector<Scalar> **x);
		void project_sub(GammaVector<Scalar> *x_sub);
		void sort_bubble(GammaVector<Scalar> *x);
		
	public:
		QPSolver(Data*, Gamma *, Theta *, Scalar);
		void init();
		void finalize();
		void compute_b();
		void solve();
		Scalar get_function_value();
		void print();
		void print(int nmb_of_spaces);
		void print(int nmb_of_spaces, bool print_A_sub);

		int get_T();
		int get_K();
		int get_dim();


};



#endif
