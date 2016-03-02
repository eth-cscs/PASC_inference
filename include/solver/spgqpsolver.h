#ifndef SPGQPSOLVER_H
#define	SPGQPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "solver/qpsolver.h"
#include "data/qpdata.h"
#include "result/qpresult.h"

#define SPGQPSOLVER_DEFAULT_MAXIT 1000;
#define SPGQPSOLVER_DEFAULT_EPS 0.0001;

#define SPGQPSOLVER_DEFAULT_M 10;
#define SPGQPSOLVER_DEFAULT_GAMMA 0.9;
#define SPGQPSOLVER_DEFAULT_SIGMA2 1.0;
#define SPGQPSOLVER_DEFAULT_ALPHAINIT 1.0;

namespace pascinference {

/* settings */
class SPGQPSolverSetting : public QPSolverSetting {
	public:
		int maxit; /* max number of iterations */
		double eps; /* precision */
		int debug_mode; /* print info about the progress */
		int m; /* size of fs */
		double gamma; 
		double sigma2;
		double alphainit; /* initial step-size */

		SPGQPSolverSetting() {
			maxit = SPGQPSOLVER_DEFAULT_MAXIT;
			eps = SPGQPSOLVER_DEFAULT_EPS;
			debug_mode = DEBUG_MODE;

			m = SPGQPSOLVER_DEFAULT_M;
			gamma = SPGQPSOLVER_DEFAULT_GAMMA;
			sigma2 = SPGQPSOLVER_DEFAULT_SIGMA2;
			alphainit = SPGQPSOLVER_DEFAULT_ALPHAINIT;

		};
		~SPGQPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  SPGQPSolverSettings:" << std::endl;
			output << "   - maxit:      " << maxit << std::endl;
			output << "   - eps:        " << eps << std::endl;
			output << "   - debug_mode: " << debug_mode << std::endl;

			output << "   - m:          " << m << std::endl;
			output << "   - gamma:      " << gamma << std::endl;
			output << "   - sigma2:     " << sigma2 << std::endl;
			output << "   - alphainit:  " << alphainit << std::endl;

		};
		
};


/* for manipulation with fs - function values for generalized Armijo condition */
class SPGQPSolver_fs {
	private:
		int m; /* the length of list */
		std::list<double> fs_list; /* the list with function values */
		
	public: 
		SPGQPSolver_fs(int new_m);
		void init(double fx);
		double get_max();		
		int get_size();
		void update(double new_fx);
		
		friend std::ostream &operator<<(std::ostream &output, SPGQPSolver_fs fs);
};



/* SPGQPSolver */ 
template<class VectorBase>
class SPGQPSolver: public QPSolver<VectorBase> {
	protected:
		const QPData<VectorBase> *data; /* data on which the solver operates */
		const QPResult<VectorBase> *result; /* here solver stores results */
	
		/* temporary vectors used during the solution process */
		void allocate_temp_vectors();
		void free_temp_vectors();
		GeneralVector<VectorBase> *g; /* gradient */
		GeneralVector<VectorBase> *d; /* projected gradient */
		GeneralVector<VectorBase> *Ad; /* A*d */
	
	public:
		SPGQPSolverSetting setting;

		SPGQPSolver();
		SPGQPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result); 
		~SPGQPSolver();


		void solve();
		void solve(SolverType type){};

		void print(std::ostream &output) const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {


/* ----- Solver ----- */
/* constructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(SPGQPSolver)CONSTRUCTOR" << std::endl;

	data = NULL;
	result = NULL;
	
	/* temp vectors */
	g = NULL;
	d = NULL;
	Ad = NULL;
}

template<class VectorBase>
SPGQPSolver<VectorBase>::SPGQPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result){
	data = &new_data;
	result = &new_result;
	
	/* allocate temp vectors */
	allocate_temp_vectors();
	
}


/* destructor */
template<class VectorBase>
SPGQPSolver<VectorBase>::~SPGQPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(SPGQPSolver)DESTRUCTOR" << std::endl;

	/* free temp vectors */
	free_temp_vectors();
}

/* prepare temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::allocate_temp_vectors(){
	GeneralVector<VectorBase> *pattern = data->b; /* I will allocate temp vectors subject to linear term */

	g = new GeneralVector<VectorBase>(*pattern);
	d = new GeneralVector<VectorBase>(*pattern);
	Ad = new GeneralVector<VectorBase>(*pattern);	
	
}

/* destroy temp_vectors */
template<class VectorBase>
void SPGQPSolver<VectorBase>::free_temp_vectors(){
	free(g);
	free(d);
	free(Ad);
	
}


/* print info about problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(SPGQPSolver)FUNCTION: print" << std::endl;

	output << " SPGQPSolver" << std::endl;
	
	/* print settings */
	output << setting;
		
}

/* solve the problem */
template<class VectorBase>
void SPGQPSolver<VectorBase>::solve() {
	if(DEBUG_MODE >= 100) std::cout << "(SPGQPSolver)FUNCTION: solve" << std::endl;

	/* I don't want to write (*x) as a vector, therefore I define following pointer types */
	typedef GeneralVector<VectorBase> (&pVector);
//	typedef GeneralMatrix<VectorBase> (&pMatrix);

	/* pointers to data */
//	pMatrix A = *(data->A);
//	pVector b = *(data->b);
//	pVector x0 = *(data->x0);

	/* pointer to result */
	pVector x = *(result->x);

	/* auxiliary vectors */
//	pVector g = *(this->g); /* gradient */
//	pVector d = *(this->d); /* A-conjugate vector */
//	pVector Ad = *(this->Ad); /* A*p */

//	x = x0; /* set approximation as initial */

	std::cout << "I am solving the problem, it will be fun!" << std::endl;

//	int it = 0; /* number of iterations */
//	int hessmult = 0; /* number of hessian multiplications */

//	double fx; /* function value */
//	SPGQPSolver_fs fs(setting.m); /* store function values for generalized Armijo condition */
//	double fx_max; /* max(fs) */
//	double xi, beta_bar, beta_hat, beta; /* for Armijo condition */
//	double dd; /* dot(d,d) */
//	double gd; /* dot(g,d) */
//	double dAd; /* dot(Ad,d) */
//	double alpha_bb; /* BB step-size */

	/* initial step-size */
//	alpha_bb = setting.alphainit;
	
	data->feasibleset->project(x);

	
}



/* ---------- SPGQPSolver_fs -------------- */

/* constructor */
SPGQPSolver_fs::SPGQPSolver_fs(int new_m){
	this->m = new_m;
}

/* init the list with function values using one initial fx */
void SPGQPSolver_fs::init(double fx){
	this->fs_list.resize(this->m, fx);
}

/* get the size of the list */
int SPGQPSolver_fs::get_size(){
	return this->m;
}

/* get the value of max value in the list */
double SPGQPSolver_fs::get_max(){
	std::list<double>::iterator it;
	it = std::max_element(this->fs_list.begin(), this->fs_list.end());
	return *it;
}

/* update the list by new value - pop the first and push the new value (FIFO) */
void SPGQPSolver_fs::update(double new_fx){
	this->fs_list.pop_back();
	this->fs_list.push_front(new_fx);
}

/* print the content of the list */
std::ostream &operator<<(std::ostream &output, SPGQPSolver_fs fs)
{
	int j, list_size;
	std::list<double>::iterator it; /* iterator through list */

	it = fs.fs_list.begin();
	list_size = fs.fs_list.size(); /* = m? */
	
	output << "[ ";
	/* for each component go throught the list */
	for(j=0;j<list_size;j++){
		output << *it;
		if(j < list_size-1){ 
				/* this is not the last node */
				output << ", ";
				it++;
		}
	}
	output << " ]";
			
	return output;
}



} /* end namespace */

#endif
