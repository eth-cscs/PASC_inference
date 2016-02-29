#ifndef DATA_H
#define	DATA_H

#include "common.h"

namespace pascinference {

class Data {
	protected:
		int dim; /* number of data components */
		int T; /* length of time serie */

	public:
		DataVector data_vec; /* data vector, TODO: should be private */

		void init(int dim, int T); /* TODO: should be constructor */
		void finalize(); /* TODO: should be destructor */

		void print();
//		void get_covtrace(PetscScalar *covtrace);
		
		virtual void generate();
		
		/* GET functions */
		int get_dim();
		int get_T();
		DataVector get_data_vec();
	
};

class Data_kmeans : public Data {
	private:
		void my_mvnrnd_D2(Scalar *mu, Scalar *covariance, Scalar *value1, Scalar *value2);
		void get_problem_value1(Scalar *value1_out, Scalar *value2_out);
		void get_problem_value2(Scalar *value1_out, Scalar *value2_out);
		void get_problem_value3(Scalar *value1_out, Scalar *value2_out);
	
	public:
		void generate();
	
};


} /* end of namespace */


#endif
