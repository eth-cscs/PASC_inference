#include "pascinference.h"

#include "matrix/localdense.h"
#include "vector/localvector.h"

#include <fstream>

namespace pascinference {
	namespace example {
		template<class VectorBase>
		class Varx2D {
			private:
				/** @brief get random number 1D 
				* 
				* return 1D random number
				* 
				* @param mu mean value 
				* @param dev standart deviation
				* @param value1 random value
				*/ 
				static void my_mvnrnd_D1(double mu, double dev, double *value1);
			
				static int get_active_cluster(int t);
				
			public:
			
				/** @brief sample generator 
				* 
				* @param T length of time-serie 
				* @param K number of clusters
				* @param mu array (of length K) with mean values for generating data 
				* @param dev array (of length K) with standart deviations for generating data
				* @param datavector output datavector with random values 
				*/ 
				static void generate(int T, int xmem, GeneralVector<VectorBase> *datavector);

		};
		
	} /* end of namespace example */
}/* end of namespace pascinference */

namespace pascinference {
	namespace example {
		template<class VectorBase>
		void Varx2D<VectorBase>::my_mvnrnd_D1(double mu, double dev, double *value1){
			double r1; 

			/* Compute normally distributed random number */
			r1 = rand()/(double)(RAND_MAX);

			/* compute output values */
			/* y = L*randn(2,1) + mean */
			*value1 = r1*dev + mu;
		}		

		template<class VectorBase>
		int Varx2D<VectorBase>::get_active_cluster(int t){
			int active_cluster = 0; 
			if((t>143 && t <=159) || (t>246 && t <= 303) || (t>346 && t <=382) || (t>433 && t <=475) || (t>577 && t <= 672) || (t>911 && t <=971) ){
				active_cluster = 1;
			} else {
				active_cluster = 0;
			} 

			return active_cluster;
		}		
		
		template<class VectorBase>
		void Varx2D<VectorBase>::generate(int T, int xmem, GeneralVector<VectorBase> *datavector) {
			/* parameters of the model */
			int K_orig = 2;
			int dimx = 2;
			
			/* mu*/
			double muK1_array[dimx] = {2.0, 6.0};
			LocalVector<VectorBase> muK1(muK1_array,dimx);
			
			double muK2_array[dimx] = {5.0, 4.0};
			LocalVector<VectorBase> muK2(muK2_array,dimx);

			/* cov */
			double covarianceK1Q1_array[dimx*dimx] = {0.03, 0.02, -0.07, 0.07};
			LocalDenseMatrix<VectorBase> covarianceK1Q1(covarianceK1Q1_array,dimx,dimx);

			double covarianceK1Q2_array[dimx*dimx] = {0.07, 0.01, -0.03, 0.06};
			LocalDenseMatrix<VectorBase> covarianceK1Q2(covarianceK1Q2_array,dimx,dimx);

			double covarianceK1Q3_array[dimx*dimx] = {-0.4, -0.2, -0.7, -0.8};
			LocalDenseMatrix<VectorBase> covarianceK1Q3(covarianceK1Q3_array,dimx,dimx);

			double covarianceK2Q1_array[dimx*dimx] = {0.03, -0.06, 0.07, 0.04};
			LocalDenseMatrix<VectorBase> covarianceK2Q1(covarianceK2Q1_array,dimx,dimx);

			double covarianceK2Q2_array[dimx*dimx] = {0.03, 0.01, 0.08, -0.09};
			LocalDenseMatrix<VectorBase> covarianceK2Q2(covarianceK2Q2_array,dimx,dimx);

			double covarianceK2Q3_array[dimx*dimx] = {0.1, -0.3, -0.2, -0.4};
			LocalDenseMatrix<VectorBase> covarianceK2Q3(covarianceK2Q3_array,dimx,dimx);

			/* first x (in mem) */
			double x0_array[dimx] = {0.3, -0.5};
			LocalVector<VectorBase> x0(x0_array,dimx);
	
			double x1_array[dimx] = {0.7, 0.1};
			LocalVector<VectorBase> x1(x1_array,dimx);

			double x2_array[dimx] = {0.1, -0.9};
			LocalVector<VectorBase> x2(x2_array,dimx);


/*	double ut_array[T+xmem];
	for(int i=0;i<T+xmem;i++){
		ut[i] = (-2.0/(double)(T+xmem)) * (i+1);
	}
*/

/*	double gamma_orig[T*K];
	for(int t=0;t<T;t++){
		if((t>143 && t <=159) || (t>246 && t <= 303) || (t>346 && t <=382) || (t>433 && t <=475) || (t>577 && t <= 672) || (t>911 && t <=971) ){
			gamma_orig[t] = 0;
			gamma_orig[t+T] = 1;
		} else {
			gamma_orig[t] = 1;
			gamma_orig[t+T] = 0;
		} 
	}
*/

		double noise_sigma = 0.0005;

		typedef GeneralVector<VectorBase> (&pVector);
		pVector data = *datavector;

		int active_cluster;
		for(int t=0;t < T;t++){ /* go through whole time series */
			active_cluster = get_active_cluster(t);
			
//			my_mvnrnd_D2(mu[k], covariance[k], &random_value1, &random_value2);
//			data(t) = random_value1;
//			data(T+t) = random_value2;
		}

	
	}
		
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
