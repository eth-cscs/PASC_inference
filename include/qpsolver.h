#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

#include "common.h"
#include "gamma.h"
#include "theta.h"
#include "data.h"


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
		

		
	public:
		QPSolver(Data*, Gamma *, Theta *, Scalar);
		void init();
		void finalize();
		void solve();
		Scalar get_function_value();
		void print();
		void print(int nmb_of_spaces);
		void print(int nmb_of_spaces, bool print_A_sub);

};



#endif
