#ifndef GAMMA_H
#define	GAMMA_H

class Gamma;
class QPSolver;

#include "common.h"
#include "data.h"
#include "theta.h"
#include "qpsolver.h"

class Gamma {
		int dim; /* number of gamma components = K */

		int global_size; /* global data size with overlap */
		int local_size; /* local data size with overlap */

		int local_begin, local_end; /* ownership range */

		PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	
		PetscRandom rnd; /* random numbers generator */
		
	public:
		Vec *gamma_vecs; /* array with data vectors, TODO: should be private */

		void init(Data, int); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print(PetscViewer v);

		void prepare_random();
		void prepare_uniform();		
		void prepare_fixed();		

		void compute(QPSolver *qp_solver, Data data, Theta theta);

		/* GET functions */
		int get_local_size();
		int get_local_begin();
		int get_local_end();
		int get_global_size();
		int get_dim();

		// TODO: this should be somewhere else, Model?
		void compute_g(Vec g, Data *data, Theta *theta);
		void compute_gk(Vec g, Data *data, Theta *theta, int k);
	
};

#endif

