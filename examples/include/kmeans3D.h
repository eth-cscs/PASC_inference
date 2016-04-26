#include "pascinference.h"
#include <fstream>


namespace pascinference {
	namespace example {
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
				template<class VectorBase>
				static void generate(int T, int K, double *mu, double *covariance, GeneralVector<VectorBase> *datavector);

				/** @brief save results into VTK file 
				* 
				* take results of the problem and save them into VTK file, which could be opened for example in ParaView
				* 
				* 
				*/
				template<class VectorBase>
				static void saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector);

				template<class VectorBase>
				static void saveVTK(std::string name_of_file, std::string extension, int T, int Knum, int *K, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector);


		};
		
	} /* end of namespace example */
}/* end of namespace pascinference */

namespace pascinference {
	namespace example {
		void KMeans3D::my_mvnrnd_D2(double *mu, double *covariance, double *value1, double *value2){
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

		void KMeans3D::my_mvnrnd_D3(double *mu, double *diag_covariance, double *value1, double *value2, double *value3){
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
		void KMeans3D::generate(int T, int K, double *mu, double *covariance, GeneralVector<VectorBase> *datavector_in) {
			int xdim = 3;
			double random_value1, random_value2, random_value3; 
			int t, k;
			int myK;

			typedef GeneralVector<VectorBase> (&pVector);
			pVector datavector = *datavector_in;

			int clusterT = ceil((double)T/(double)K);
			int Tlocal = ((double)datavector.local_size())/((double)xdim);
			int t_begin, t_end; /* ownership */ 

			TRY( VecGetOwnershipRange(datavector.get_vector(), &t_begin, &t_end) );
			t_begin = ((double)t_begin)/((double)xdim);
			t_end = ((double)t_end)/((double)xdim);
	
			double *muK, *covarianceK;
			double *datavector_arr;
									
			/* generate local random data */
			TRY( VecGetArray(datavector.get_vector(),&datavector_arr) );
			for(t=0;t < Tlocal;t++){
				/* which cluster am I computing? */
				for(k=0;k<K;k++){
					if( ((t_begin+t) >= k*clusterT) && ((t_begin+t) < (k+1)*clusterT) ){
						myK = k;
					}
				}
				
				muK = &mu[myK*xdim]; /* get subarray */
				covarianceK = &covariance[myK*xdim]; /* get subarray */
				
				my_mvnrnd_D3(muK, covarianceK, &random_value1, &random_value2, &random_value3);
				
				datavector_arr[xdim*t] = random_value1;
				datavector_arr[xdim*t+1] = random_value2;
				datavector_arr[xdim*t+2] = random_value3;
				
			}
			TRY( VecRestoreArray(datavector.get_vector(),&datavector_arr) );

		}
		
		//TODO: this should be written for every type of vectorbase
		template<class VectorBase>
		void KMeans3D::saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector_in, GeneralVector<VectorBase> *gammavector_in){
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



		template<>
		void KMeans3D::saveVTK(std::string name_of_file, std::string extension, int T, int Knum, int *Karr, GeneralVector<PetscVector> *datavector_in, GeneralVector<PetscVector> *gammavector_in){
			Timer timer_saveVTK; 
			timer_saveVTK.restart();
			timer_saveVTK.start();

			/* from input pointers to classical vectors (just simplification of notation) */
			typedef GeneralVector<PetscVector> (&pVector);
			pVector datavector = *datavector_in;
			pVector gammavector = *gammavector_in;
	
			int t,k;
			int my_rank = GlobalManager.get_rank();
			int K = Karr[my_rank];
			int Tlocal = datavector.local_size()/3.0; /* xdim = 3 */
			
			/* filename */
			std::ostringstream oss_name_of_file;
			
			/* to manipulate with file */
			std::ofstream myfile;
	
			/* write to the name of file */
			oss_name_of_file << name_of_file << "_" << my_rank << extension;

			/* open file to write */
			myfile.open(oss_name_of_file.str().c_str());

			/* write header to file */
			myfile << "# vtk DataFile Version 3.1\n";
			myfile << "PASCInference: Kmeans2D solution, K=" << K << "\n";
			myfile << "ASCII\n";
			myfile << "DATASET UNSTRUCTURED_GRID\n";

			/* points - coordinates */
			myfile << "POINTS " << T << " FLOAT\n";

			/* sequential part follows, therefore close the file for now */
			myfile.close();
			TRY(PetscBarrier(NULL));
			
			/* each processor will write into one folder */ //TODO: this is slow :(
			std::ostringstream oss_name_of_other_file;
			std::ofstream otherfile;
			double *datavector_arr;
			for(int idproc = 0; idproc < GlobalManager.get_size(); idproc++){
				if(idproc == my_rank){
					TRY( VecGetArray(datavector.get_vector(),&datavector_arr) )

					
					/* now I will write into files */
					for(int idproc2 = 0; idproc2 < GlobalManager.get_size(); idproc2++){
						oss_name_of_other_file << name_of_file << "_" << idproc2 << extension;
						otherfile.open(oss_name_of_other_file.str().c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
						
						for(t=0;t < Tlocal;t++){
							otherfile << datavector_arr[3*t] << " "; /* x */
							otherfile << datavector_arr[3*t+1] << " "; /* y */
							otherfile << datavector_arr[3*t+2] << "\n"; /* z */
						}
						
						otherfile.close();
						oss_name_of_other_file.str("");
						oss_name_of_other_file.clear();
					}

					TRY( VecRestoreArray(datavector.get_vector(),&datavector_arr) )
				}

				TRY(PetscBarrier(NULL));
			}

			/* parallel part continues */
			myfile.open(oss_name_of_file.str().c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
			myfile << "\n";

			/* values in points */
			myfile << "POINT_DATA " << T << "\n";
			/* prepare vector with idx of max values */

			/* get local vector of gamma */
			double *gammavector_arr;
			TRY( VecGetArray(gammavector.get_vector(),&gammavector_arr) );

			double gamma_max[gammavector.local_size()];
			int gamma_max_idx[gammavector.local_size()];
			double temp;
			/* initialization of arrays */
			for(t=0;t<T;t++){
				gamma_max[t] = 0.0;
				gamma_max_idx[t] = 0;
			}
	
			for(k=0;k<K;k++){
				/* write gamma_k */
				myfile << "SCALARS gamma_" << k << " float 1\n";
				myfile << "LOOKUP_TABLE default\n";
				for(t=0;t<T;t++){
					myfile << gammavector_arr[k*T + t] << "\n";

					/* update maximum */
					if(gammavector_arr[k*T + t] > gamma_max[t]){
						temp = gammavector_arr[k*T+t];
						gamma_max[t] = temp;
						gamma_max_idx[t] = k;
					}
				}
			}

			TRY( VecRestoreArray(gammavector.get_vector(),&gammavector_arr) )

			/* store gamma values */
			myfile << "SCALARS gamma_max_id float 1\n";
			myfile << "LOOKUP_TABLE default\n";
			for(t=0;t<T;t++){
				myfile << gamma_max_idx[t] << "\n";
			}

			/* store proc_id */
			myfile << "SCALARS proc_id float 1\n";
			myfile << "LOOKUP_TABLE default\n";
			/* sequential part follows, therefore close the file for now */
			myfile.close();
			TRY(PetscBarrier(NULL));
			
			/* each processor will write into one folder his id */ //TODO: this is slow :(
			for(int idproc = 0; idproc < GlobalManager.get_size(); idproc++){
				if(idproc == my_rank){
					/* now I will write into files */
					for(int idproc2 = 0; idproc2 < GlobalManager.get_size(); idproc2++){
						oss_name_of_other_file << name_of_file << "_" << idproc2 << extension;
						otherfile.open(oss_name_of_other_file.str().c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
						
						for(t=0;t < Tlocal;t++){
							otherfile << my_rank << "\n"; /* x */
						}
						
						otherfile.close();
						oss_name_of_other_file.str("");
						oss_name_of_other_file.clear();
					}
				}

				TRY(PetscBarrier(NULL));
			}

			timer_saveVTK.stop();
			coutAll <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;

		}		
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
