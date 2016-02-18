#ifndef QPPROBLEM_H
#define	QPPROBLEM_H

class QPSolver;

#include "common.h"

#include "operations.h"
#include "projection.h"
#include <list>

/* for manipulation with fs - function values for generalized Armijo condition */
class QPSolver_fs {
	private:
		int m; /* the length of list */
		std::list<Scalar> fs_list; /* the list with function values */
		
	public: 
		QPSolver_fs(int new_m);
		void init(Scalar fx);
		Scalar get_max();		
		int get_size();
		void update(Scalar new_fx);
		
		friend std::ostream &operator<<(std::ostream &output, QPSolver_fs fs);
	
};

class QPSolver {
	protected:
		int T;
		int K;

		Scalar eps_sqr;

		int it;
		int it_all;
		int hessmult;
		int hessmult_all;
		
		/* data of QP problem */
		GammaVector g; /* gradient */
		GammaVector d; /* projected gradient */
		GammaVector Ad; /* A*ds */
		
		Timer timer_projection; /* the sum of time necessary to perform projections */
		Timer timer_matmult; /* the sum of time necessary to perform matrix multiplication */
		Timer timer_dot; /* the sum of time necessary to compute dot_products */
		Timer timer_update; /* total time of vector updates */
		Timer timer_stepsize; /* total time of step-size computation */
		Timer timer_fs; /* total time of manipulation with fs vector during iterations */
		Timer timer_total; /* total time of SPG algorithm */
		
	public:
		GammaVector b; /* rhs */ // TODO: this should be private

		void init(int T, int K, Scalar eps_sqr);
		void finalize();

		void solve(GammaVector &x);
		Scalar get_function_value(GammaVector &x);
		Scalar get_function_value(GammaVector &x, bool use_gradient);

		void print();
		void print_timers();
		void print(int nmb_of_spaces);

		int get_T();
		int get_K();

		int get_it();
		int get_it_all();
		int get_hessmult();
		int get_hessmult_all();
		double get_time_projection();
		double get_time_matmult();
		double get_time_dot();
		double get_time_update();
		double get_time_stepsize();
		double get_time_fs();
		double get_time_total();

};



#endif
