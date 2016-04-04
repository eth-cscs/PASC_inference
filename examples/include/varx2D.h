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
				static void generate(int T, int xmem, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *u);

				/** @brief save results into VTK file 
				* 
				* take results of the problem and save them into VTK file, which could be opened for example in ParaView
				* 
				*/
				static void saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector);


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
		void Varx2D<VectorBase>::generate(int T, int xmem, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *uvector) {
			/* parameters of the model */
			int dimx = 2;
			
			/* mu*/
			double muK1_array[dimx] = {2.0, 6.0};
			double muK2_array[dimx] = {5.0, 4.0};

			/* A in column major order */
			double A_K1Q1_array[dimx*dimx] = {0.03, 0.02, -0.07, 0.07};
			double A_K1Q2_array[dimx*dimx] = {0.07, 0.01, -0.03, 0.06};
			double A_K1Q3_array[dimx*dimx] = {-0.4, -0.2, -0.7, -0.8};

			double A_K2Q1_array[dimx*dimx] = {0.03, -0.06, 0.07, 0.04};
			double A_K2Q2_array[dimx*dimx] = {0.03, 0.01, 0.08, -0.09};
			double A_K2Q3_array[dimx*dimx] = {0.1, -0.3, -0.2, -0.4};

			/* B */
			double B_K1Q0_array[dimx*dimx] = {0.4, 0.1};
			double B_K1Q1_array[dimx*dimx] = {-0.3, 0.6};

			double B_K2Q0_array[dimx*dimx] = {0.4, 0.1};
			double B_K2Q1_array[dimx*dimx] = {-0.1, 0.2};

			double noise_sigma = 0.0005;

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
				active_cluster = get_active_cluster(t);

				/* create time series by recursive formula data(t) = A_Q1*data(t-1) + A_Q2*data(t-2) + A_Q3*data(t-3) + B_Q0*u(t) + B_Q1*u(t-1) + noise; */
				if(active_cluster == 0){
					data(t         ) = 	muK1_array[0]
										+ (A_K1Q1_array[0]*data(t-1)+A_K1Q1_array[2]*data(t-1+datasize))
										+ (A_K1Q2_array[0]*data(t-2)+A_K1Q2_array[2]*data(t-2+datasize))
										+ (A_K1Q3_array[0]*data(t-3)+A_K1Q3_array[2]*data(t-3+datasize))
										+ (B_K1Q0_array[0]*u(t) + B_K1Q1_array[0]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
										
					data(t+datasize) = 	muK1_array[1]
										+ (A_K1Q1_array[1]*data(t-1)+A_K1Q1_array[3]*data(t-1+datasize))
										+ (A_K1Q2_array[1]*data(t-2)+A_K1Q2_array[3]*data(t-2+datasize))
										+ (A_K1Q3_array[1]*data(t-3)+A_K1Q3_array[3]*data(t-3+datasize))
										+ (B_K1Q0_array[1]*u(t) + B_K1Q1_array[1]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
					
				}
			
				if(active_cluster == 1){
					data(t         ) = 	muK2_array[0]
										+ (A_K2Q1_array[0]*data(t-1)+A_K2Q1_array[2]*data(t-1+datasize))
										+ (A_K2Q2_array[0]*data(t-2)+A_K2Q2_array[2]*data(t-2+datasize))
										+ (A_K2Q3_array[0]*data(t-3)+A_K2Q3_array[2]*data(t-3+datasize))
										+ (B_K2Q0_array[0]*u(t) + B_K2Q1_array[0]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
										
					data(t+datasize) = 	muK1_array[1]
										+ (A_K2Q1_array[1]*data(t-1)+A_K2Q1_array[3]*data(t-1+datasize))
										+ (A_K2Q2_array[1]*data(t-2)+A_K2Q2_array[3]*data(t-2+datasize))
										+ (A_K2Q3_array[1]*data(t-3)+A_K2Q3_array[3]*data(t-3+datasize))
										+ (B_K2Q0_array[1]*u(t) + B_K2Q1_array[1]*u(t-1))
										+ noise_sigma*rand()/(double)(RAND_MAX);
					
				}


			}
	
		}


		//TODO: this should be written for every type of vectorbase
		template<class VectorBase>
		void Varx2D<VectorBase>::saveVTK(std::string name_of_file, int T, int K, GeneralVector<VectorBase> *datavector_in, GeneralVector<VectorBase> *gammavector_in){
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
					myfile << datavector.get(t) << " "; /* x */
					myfile << datavector.get(T+t) << " "; /* y */
					myfile << " 0.0\n"; /* z */
			}
			myfile << "\n";

			//TODO: store gamma right here

			/* close file */
			myfile.close();

			timer_saveVTK.stop();
			coutMaster <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;

		}


		
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
