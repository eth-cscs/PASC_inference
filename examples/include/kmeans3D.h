#include "pascinference.h"
#include <fstream>


namespace pascinference {
	namespace example {
		template<class VectorBase>
		class KMeans3D {
			private:
				/** @brief get random number 2D 
				* 
				* return 2D random number
				* 
				* @param mu array (of length 2) with mean value 
				* @param diagonal covariance array (of length 2) with diagonal covariance matrix
				* @param value1 x-coordinate of random point 
				* @param value2 y-coordinate of random point 
				*/ 
				static void my_mvnrnd_D2(double *mu, double *covariance, double *value1, double *value2);
			
				/** @brief get random number 3D 
				* 
				* return 3D random number
				* 
				* @param mu array (of length 3) with mean value 
				* @param diagonal covariance array (of length 3) with diagonal covariance matrix
				* @param value1 x-coordinate of random point 
				* @param value2 y-coordinate of random point 
				* @param value3 z-coordinate of random point 
				*/ 
				static void my_mvnrnd_D3(double *mu, double *diag_covariance, double *value1, double *value2, double *value3);
			
			public:
			
				/** @brief sample generator 
				* 
				* @param T length of time-serie 
				* @param K number of clusters (i.e. solution)
				* @param mu array (of length K*dimx) with mean values for generating data 
				* @param covariance array (of length K*dimx) of matrices of size dimx*dimx with covariance matrices for generating data
				* @param datavector output datavector with random values 
				*/ 
				static void generate(int T, int K, double *mu, double *covariance, GeneralVector<VectorBase> *datavector);

				/** @brief save results into VTK file 
				* 
				* take results of the problem and save them into VTK file, which could be opened for example in ParaView
				* 
				* 
				*/
				static void saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector);

		};
		
	} /* end of namespace example */
}/* end of namespace pascinference */

namespace pascinference {
	namespace example {
		template<class VectorBase>
		void KMeans3D<VectorBase>::my_mvnrnd_D2(double *mu, double *covariance, double *value1, double *value2){
			double L[4];
			double r1, r2, r1n, r2n; 
	
			double R, c, s;

			/* Compute normally distributed random numbers via Box-Muller */
			r1 = rand()/(double)(RAND_MAX);
			r1 = 1.-r1; /* to change from [0,1) to (0,1], which we need for the log */

			r2 = rand()/(double)(RAND_MAX);
			R = sqrt(-2.*log(r1));
			c = cos(2.*M_PI*r2);
			s = sin(2.*M_PI*r2);

			/* compute normal distributed random values */
			r1n = R*c;
			r2n = R*s;
	
			/* choleski decomposition of SPD covariance matrix */
			L[0] = sqrt(covariance[0]);
			L[1] = 0;
			L[2] = covariance[1]/L[0];
			L[3] = sqrt(covariance[3] - L[2]*L[2]);

			/* compute output values */
			/* y = L*randn(2,1) + mean */
			*value1 = L[0]*r1n + L[1]*r2n + mu[0];
			*value2 = L[2]*r1n + L[3]*r2n + mu[1];

		}		

		template<class VectorBase>
		void KMeans3D<VectorBase>::my_mvnrnd_D3(double *mu, double *diag_covariance, double *value1, double *value2, double *value3){
			double r1, r2, r3, r4;

			double mu12[2] = {mu[0],mu[1]};
			double mu34[2] = {mu[2],0.0};
			
			double diag_covariance12[4] = {diag_covariance[0],0.0,0.0,diag_covariance[1]};
			double diag_covariance34[4] = {diag_covariance[2],0.0,0.0,1.0};

			my_mvnrnd_D2(mu12, diag_covariance12, &r1, &r2);
			my_mvnrnd_D2(mu34, diag_covariance34, &r3, &r4);

			/* compute output values */
			*value1 = r1;
			*value2 = r2;
			*value3 = r3;

			/* discart r4 */
		}		
		
		template<class VectorBase>
		void KMeans3D<VectorBase>::generate(int T, int K, double *mu, double *covariance, GeneralVector<VectorBase> *datavector) {
			int dimx = 3;
			double random_value1, random_value2, random_value3; 
			int t, k;

			typedef GeneralVector<VectorBase> (&pVector);
			pVector data = *datavector;

			int clusterT = ceil((double)T/(double)K);
	
			double *muK, *covarianceK;
	
			/* generate random data */
			for(k=0;k<K;k++){ /* go through clusters */ 
				muK = &mu[k*dimx]; /* get subarray */
				covarianceK = &covariance[k*dimx]; /* get subarray */
				for(t=k*clusterT;t < std::min((k+1)*clusterT,T);t++){ /* go through part of time serie corresponding to this cluster */
					my_mvnrnd_D3(muK, covarianceK, &random_value1, &random_value2, &random_value3);
		
					data(t) = random_value1;
					data(T+t) = random_value2;
					data(2*T+t) = random_value3;
				}
			}
		}
		
		//TODO: this should be written for every type of vectorbase
		template<class VectorBase>
		void KMeans3D<VectorBase>::saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector_in, GeneralVector<VectorBase> *gammavector_in){
			Timer timer_saveVTK; 
			timer_saveVTK.restart();
			timer_saveVTK.start();

			/* from input pointers to classical vectors (just simplification of notation) */
			typedef GeneralVector<VectorBase> (&pVector);
			pVector datavector = *datavector_in;
			pVector gammavector = *gammavector_in;
	
			int t,k;
			
			/* filename */
			std::ostringstream oss_name_of_file;
			
			/* to manipulate with file */
			std::ofstream myfile;
	
			/* write to the name of file */
			oss_name_of_file << name_of_file;

			/* open file to write */
			myfile.open(oss_name_of_file.str().c_str());

			/* write header to file */
			myfile << "# vtk DataFile Version 3.1\n";
			myfile << "PASCInference: Kmeans2D solution\n";
			myfile << "ASCII\n";
			myfile << "DATASET UNSTRUCTURED_GRID\n";

			/* points - coordinates */
			myfile << "POINTS " << T << " FLOAT\n";
			for(t=0;t < T;t++){
					myfile << datavector.get(t) << " "; /* x */
					myfile << datavector.get(T+t) << " "; /* y */
					myfile << datavector.get(2*T+t) << "\n"; /* z */
			}
			myfile << "\n";

			/* values in points */
			myfile << "POINT_DATA " <<  T << "\n";
			/* prepare vector with idx of max values */

			GeneralVector<VectorBase> gamma_max(gammavector);
			double temp;
			gamma_max(gall) = 0.0;
			GeneralVector<VectorBase> gamma_max_idx(gammavector); // TODO: use general host vector
	
			gamma_max_idx(gall) = 0;
			for(k=0;k<K;k++){
				/* write gamma_k */
				myfile << "SCALARS gamma_" << k << " float 1\n";
				myfile << "LOOKUP_TABLE default\n";
				for(t=0;t<T;t++){
					myfile << gammavector.get(k*T + t) << "\n";

					/* update maximum */
					if(gammavector(k*T+t) > gamma_max(t)){
						temp = gammavector.get(k*T+t);
						gamma_max(t) = temp;
				
//						gamma_vec(k*T+t) = gamma_vec(gamma_max_idx(t)*T+t);
						gamma_max_idx(t) = k;
					}
				}
			}


			/* store gamma values */
			myfile << "SCALARS gamma_max_id float 1\n";
			myfile << "LOOKUP_TABLE default\n";
			for(t=0;t<T;t++){
				myfile << (int)gamma_max_idx.get(t) << "\n";
			}

			/* close file */
			myfile.close();

			timer_saveVTK.stop();
			coutMaster <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;

		}
		
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
