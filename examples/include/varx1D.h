#include "pascinference.h"
#include <fstream>

namespace pascinference {
	namespace example {
		template<class VectorBase>
		class Varx1D {
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
			
			public:
			
				/** @brief sample generator 
				* 
				* @param T length of time-serie 
				* @param K number of clusters
				* @param mu array (of length K) with mean values for generating data 
				* @param dev array (of length K) with standart deviations for generating data
				* @param datavector output datavector with random values 
				*/ 
				static void generate(int T, int K, double *mu, double *dev, GeneralVector<VectorBase> *datavector);

		};
		
	} /* end of namespace example */
}/* end of namespace pascinference */

namespace pascinference {
	namespace example {
		template<class VectorBase>
		void Varx1D<VectorBase>::my_mvnrnd_D1(double mu, double dev, double *value1){
			double r1; 

			/* Compute normally distributed random number */
			r1 = rand()/(double)(RAND_MAX);

			/* compute output values */
			/* y = L*randn(2,1) + mean */
			*value1 = r1*dev + mu;
		}		
		
		template<class VectorBase>
		void Varx1D<VectorBase>::generate(int T, int K, double *mu, double *dev, GeneralVector<VectorBase> *datavector) {
			int t;
			int k;
			double random_value1; 

			typedef GeneralVector<VectorBase> (&pVector);
			pVector data = *datavector;

			int clusterT = ceil((double)T/(double)K); /* length of the cluster */
	
			/* generate random data */
			for(k=0;k<K;k++){ /* go through clusters */ 
				for(t=k*clusterT;t < std::min((k+1)*clusterT,T);t++){ /* go through part of time serie corresponding to this cluster */
					my_mvnrnd_D1(mu[k], dev[k], &random_value1);
		
					data(t) = random_value1;
				}
			}
		}
		
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
