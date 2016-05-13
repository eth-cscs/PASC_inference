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
				* @param xmem 
				* @param datavector output datavector with random values 
				* @param u 
				*/ 
				static void generate(int T, int xmem, double noise_sigma, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *u);

				static void generate(int T, int xmem, double noise_sigma, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *u, GeneralVector<VectorBase> *thetavector_solution, GeneralVector<VectorBase> *gammavector_solution);

				/** @brief save results into VTK file 
				* 
				* take results of the problem and save them into VTK file, which could be opened for example in ParaView
				* 
				*/
				static void saveVTK(std::string name_of_file, int T, int xmem, int K, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector);


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
		void Varx2D<VectorBase>::generate(int T, int xmem, double noise_sigma, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *uvector) {
			/* parameters of the model */
			int dimx = 2;
			int K = 2;
			int umem = 1;
			
			/* theta_solution */
			int theta_size_n = 1 + dimx*xmem + (umem + 1);
			int theta_size_k = dimx * theta_size_n; 
			int theta_size = K*theta_size_k; 
			
			double theta[theta_size] =
					{ 2.0,   0.03,-0.07,  0.07,-0.03,  -0.4,-0.7,  0.4,	 -0.3,		/* K=1, n=1:  mu,A1,A2,A3,B1,B2 */
					  6.0,   0.02, 0.07,  0.01, 0.06,  -0.2,-0.8,  0.1,   0.6, 		/* K=1, n=2 */
					  5.0,   0.03, 0.07,  0.03, 0.08,   0.1,-0.2,  0.4,  -0.1,		/* K=2, n=1 */
					  4.0,  -0.06, 0.04,  0.01,-0.09,  -0.3,-0.4,  0.1,   0.2		/* K=2, n=2 */
					};

			typedef GeneralVector<VectorBase> (&pVector);
			pVector data = *datavector;
			pVector u = *uvector;

			/* first x (in mem) */
			data(0         ) = 0.3;
			data(0+(T+xmem)) = -0.5;

			data(1         ) = 0.7;
			data(1+(T+xmem)) = 0.1;

			data(2         ) = 0.1;
			data(2+(T+xmem)) = -0.9;

			/* fill all other data */
			int datasize = T+xmem;
			/* fill u */
			for(int t=0;t < T+xmem;t++){ /* go through whole time series */
				/* uvector is given, then generate some data also to this vector */
				u(t) = (-2.0/(double)(T+xmem)) * (t+1);
			}

			int active_cluster; /* number of active cluster */
			for(int t=xmem;t < T+xmem;t++){ /* go through whole time series */
				active_cluster = get_active_cluster(t-xmem);

				/* create time series by recursive formula data(t) = A_Q1*data(t-1) + A_Q2*data(t-2) + A_Q3*data(t-3) + B_Q0*u(t) + B_Q1*u(t-1) + noise; */
				if(active_cluster == 0){
					data(t         ) = 	theta[0]
										+ (theta[1]*data(t-1) + theta[2]*data(t-1+datasize))
										+ (theta[3]*data(t-2) + theta[4]*data(t-2+datasize))
										+ (theta[5]*data(t-3) + theta[6]*data(t-3+datasize))
										+ (theta[7]*u(t) + theta[8]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
										
					data(t+datasize) = 	theta[9]
										+ (theta[10]*data(t-1) + theta[11]*data(t-1+datasize))
										+ (theta[12]*data(t-2) + theta[13]*data(t-2+datasize))
										+ (theta[14]*data(t-3) + theta[15]*data(t-3+datasize))
										+ (theta[16]*u(t) + theta[17]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
					
				}
			
				if(active_cluster == 1){
					data(t         ) = 	theta[18]
										+ (theta[19]*data(t-1) + theta[20]*data(t-1+datasize))
										+ (theta[21]*data(t-2) + theta[22]*data(t-2+datasize))
										+ (theta[23]*data(t-3) + theta[24]*data(t-3+datasize))
										+ (theta[25]*u(t) + theta[26]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
										
					data(t+datasize) = 	theta[27]
										+ (theta[28]*data(t-1) + theta[29]*data(t-1+datasize))
										+ (theta[30]*data(t-2) + theta[31]*data(t-2+datasize))
										+ (theta[32]*data(t-3) + theta[33]*data(t-3+datasize))
										+ (theta[34]*u(t) + theta[35]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
					
				}


			}
	
		}

		template<class VectorBase>
		void Varx2D<VectorBase>::generate(int T, int xmem, double noise_sigma, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *uvector, GeneralVector<VectorBase> *thetavector_solution, GeneralVector<VectorBase> *gammavector_solution) {
			/* generate data using normal way */
			generate(T, xmem, noise_sigma, datavector, uvector);

			/* fill gamma solution vector */
			int i, active_cluster;
			
			for(i=0;i<T;i++){
				active_cluster = get_active_cluster(i);
				if(active_cluster == 0){
					(*gammavector_solution)(i) = 1.0;
				}
				if(active_cluster == 1){
					(*gammavector_solution)(i+T) = 1.0;
				}
			}

			/* fill thetavector */
			int dimx = 2;
			int umem = 1;
			int K = 2;
			int theta_size_n = 1 + dimx*xmem + (umem + 1);
			int theta_size_k = dimx * theta_size_n; 
			int theta_size = K*theta_size_k; 			
			double theta[theta_size] =
					{ 2.0,   0.03,-0.07,  0.07,-0.03,  -0.4,-0.7,  0.4,	 -0.3,		/* K=1, n=1:  mu,A1,A2,A3,B1,B2 */
					  6.0,   0.02, 0.07,  0.01, 0.06,  -0.2,-0.8,  0.1,   0.6, 		/* K=1, n=2 */
					  5.0,   0.03, 0.07,  0.03, 0.08,   0.1,-0.2,  0.4,  -0.1,		/* K=2, n=1 */
					  4.0,  -0.06, 0.04,  0.01,-0.09,  -0.3,-0.4,  0.1,   0.2		/* K=2, n=2 */
					};
			for(i=0;i<theta_size;i++){
				(*thetavector_solution)(i) = theta[i];
			}		

		}

		//TODO: this should be written for every type of vectorbase
		template<class VectorBase>
		void Varx2D<VectorBase>::saveVTK(std::string name_of_file, int T, int xmem, int K, GeneralVector<VectorBase> *datavector_in, GeneralVector<VectorBase> *gammavector_in){
			Timer timer_saveVTK; 
			timer_saveVTK.restart();
			timer_saveVTK.start();

			/* from input pointers to classical vectors (just simplification of notation) */
			typedef GeneralVector<VectorBase> (&pVector);
			pVector datavector = *datavector_in;
			pVector gammavector = *gammavector_in;
	
			int t;
			
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
			myfile << "PASCInference: Varx2D solution\n";
			myfile << "ASCII\n";
			myfile << "DATASET UNSTRUCTURED_GRID\n";

			/* points - coordinates */
			myfile << "POINTS " << T << " FLOAT\n";
			for(t=0;t < T;t++){
					myfile << datavector.get(t+xmem) << " "; /* x */
					myfile << datavector.get(2*xmem+T+t) << " "; /* y */
					myfile << " 0.0\n"; /* z */
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
			int k;
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
