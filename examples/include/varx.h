#include "pascinference.h"
#include <fstream>


namespace pascinference {
	namespace example {
		class VarX {
			private:

				static void compute_next_step(double *data_out, int t_data_out, double *data_in, int t_data_in, int xdim, int xmem, double *theta, int k);
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
				static void generate(int T, int xdim, int K, int xmem, double *theta, double *xstart, int (*get_cluster_id)(int, int), GeneralVector<VectorBase> *datavector, double noise, bool scale_or_not);

				/** @brief save results into CSV file 
				* 
				*/
				template<class VectorBase>
				static void saveCSV(std::string name_of_file, int T, int xdim, int *xmem_arr, int *K_arr, GeneralVector<VectorBase> *datavector, GeneralVector<VectorBase> *gammavector, GeneralVector<VectorBase> *thetavector);

				template<class VectorBase>
				static void set_solution_gamma(int T, int xdim, int K, int xmem, int (*get_cluster_id)(int, int), GeneralVector<VectorBase> *gammavector);

				template<class VectorBase>
				static void set_solution_theta(int T, int xdim, int K, int xmem, double *theta, GeneralVector<VectorBase> *thetavector);

		};
		
	} /* end of namespace example */
}/* end of namespace pascinference */

namespace pascinference {
	namespace example {

		void VarX::compute_next_step(double *data_out, int t_data_out, double *data_in, int t_data_in, int xdim, int xmem, double *theta, int k){
			int theta_length_n = 1+xdim*xmem; /* the size of theta for one dimension, one cluster */
			int theta_start = k*theta_length_n*xdim; /* where in theta start actual coefficients */

			double Ax;

			int t_mem,n,i;
			for(n = 0; n < xdim; n++){
				/* mu */
				data_out[t_data_out*xdim+n] = theta[theta_start + n*theta_length_n]; 
				
				/* A */
				for(t_mem = 1; t_mem <= xmem; t_mem++){
					/* add multiplication with A_{t_mem} */
					Ax = 0;
					for(i = 0; i < xdim; i++){
//						coutMaster << "A" << t_mem << "_" << n << "," << i << " = " << theta[theta_start + n*theta_length_n + 1 + (t_mem-1)*xdim + i] << std::endl;
						Ax += theta[theta_start + n*theta_length_n + 1 + (t_mem-1)*xdim + i]*data_in[(t_data_in-t_mem)*xdim+i]; 
					}
					data_out[t_data_out*xdim+n] += Ax;
		
				}
			
			}
		}


#ifdef USE_PETSCVECTOR		
		template<>
		void VarX::set_solution_gamma(int T, int xdim, int K, int xmem, int (*get_cluster_id)(int, int), GeneralVector<PetscVector> *gammavector){
			LOG_FUNC_STATIC_BEGIN
			
			double *gammavector_arr;
			int k,t;
			
			TRY( VecSet(gammavector->get_vector(),0.0) );
			TRY( VecGetArray(gammavector->get_vector(),&gammavector_arr) );

			for(t=xmem;t<T;t++){
				k = (*get_cluster_id)(t, T);
				gammavector_arr[k*(T-xmem) + t - xmem] = 1;
			}

			TRY( VecRestoreArray(gammavector->get_vector(),&gammavector_arr) );
			
			LOG_FUNC_STATIC_END
		}		

		template<>
		void VarX::set_solution_theta(int T, int xdim, int K, int xmem, double *theta, GeneralVector<PetscVector> *thetavector){
			LOG_FUNC_STATIC_BEGIN

			double *thetavector_arr;
			int i;
			
			TRY( VecSet(thetavector->get_vector(),0.0) );
			TRY( VecGetArray(thetavector->get_vector(),&thetavector_arr) );

			for(i=0;i<K*xdim*(1+xdim*xmem);i++){
				thetavector_arr[i] = theta[i];
			}

			TRY( VecRestoreArray(thetavector->get_vector(),&thetavector_arr) );
			
			LOG_FUNC_STATIC_END
		}		

		template<>
		void VarX::generate(int T, int xdim, int K, int xmem, double *theta, double *xstart, int (*get_cluster_id)(int, int), GeneralVector<PetscVector> *datavector_in, double noise, bool scale_or_not) {
			LOG_FUNC_STATIC_BEGIN
			
			/* size of input */
			int theta_size = K*xdim*(1 + xdim*xmem); 
			int xstart_size = xdim*xmem;

			int nproc = GlobalManager.get_size();
			int my_rank = GlobalManager.get_rank();

			/* get ownership range */
			int t_begin, t_end, t_length; 
			TRY( VecGetOwnershipRange(datavector_in->get_vector(), &t_begin, &t_end) );
			t_begin = ((double)t_begin)/((double)xdim);
			t_end = ((double)t_end)/((double)xdim);			
			t_length = t_end - t_begin;
			
			/* scattering the tails */
			VecScatter ctx; /* for scattering xn_global to xn_local */
			IS scatter_is;
			IS scatter_is_to;
			TRY( ISCreateStride(PETSC_COMM_SELF, xstart_size, 0, 1, &scatter_is_to) );
			
			Vec xtail_global;
				TRY( VecCreate(PETSC_COMM_WORLD,&xtail_global) );
				TRY( VecSetSizes(xtail_global,xstart_size, PETSC_DECIDE) );
				TRY( VecSetFromOptions(xtail_global) );
			Vec xtail_local;
				TRY( VecCreateSeq(PETSC_COMM_SELF, xstart_size, &xtail_local) );
			double *xtail_arr;

			/* get local data array */
			double *x_arr;
			TRY( VecGetArray(datavector_in->get_vector(),&x_arr) );

			/* I suppose that t_length >= xmem */
			int rank, t, n, k;
			for(rank = 0; rank < nproc; rank++){ /* through processors - this is sequential */
				/* xstart */
				if(rank == my_rank){
					/* the start */
					for(t = 0; t < xmem; t++){
						if(rank == 0){
							/* given start */
							for(n = 0; n < xdim; n++){
								x_arr[t*xdim+n] = xstart[t*xdim+n];
							}	
						} else {
							/* compute start from obtained tail from previous rank */
							TRY( VecGetArray(xtail_local,&xtail_arr) );
							for(t = 0; t < xstart_size;t++){
								x_arr[t] = xtail_arr[t];
							}
							TRY( VecRestoreArray(xtail_global,&xtail_arr) );
						}
					}
					
					for(t = xmem; t < t_length; t++){ /* through local time */
						k = (*get_cluster_id)(t_begin + t, T);
						compute_next_step(x_arr, t, x_arr, t, xdim, xmem, theta, k);

					}

				}

				/* fill the tail vector */
				TRY( VecGetArray(xtail_global,&xtail_arr) );
				for(t = 0; t < xstart_size;t++){
					xtail_arr[t] = x_arr[(t_length-xmem)*xdim+t];
				}
				TRY( VecRestoreArray(xtail_global,&xtail_arr) );

				/* now tails are stored in global vector, scatter them to local vector */
				TRY( ISCreateStride(PETSC_COMM_WORLD, xstart_size, rank*xstart_size, 1, &scatter_is) );

				/* scatter xtail_global to xtail_local */
				TRY( VecScatterCreate(xtail_global,scatter_is,xtail_local,scatter_is_to,&ctx) );
				TRY( VecScatterBegin(ctx,xtail_global,xtail_local,INSERT_VALUES,SCATTER_FORWARD) );
				TRY( VecScatterEnd(ctx,xtail_global,xtail_local,INSERT_VALUES,SCATTER_FORWARD) );

				TRY( ISDestroy(&scatter_is) );
				TRY( VecScatterDestroy(&ctx) );

				TRY( PetscBarrier(NULL) );
			}

			/* add noise to generated data */
			for(t = 0; t < t_length; t++){
				for(n = 0; n <xdim; n++){
					x_arr[t*xdim+n] += noise*rand()/(double)(RAND_MAX);			
				}
			}


			/* restore local data array */
			TRY( VecRestoreArray(datavector_in->get_vector(),&x_arr) );

			TRY( VecDestroy(&xtail_global) );
			TRY( VecDestroy(&xtail_local) );


			/* scale data */
			if(scale_or_not){
				double max_value;
				TRY( VecMax(datavector_in->get_vector(), NULL, &max_value) );
				TRY( VecScale(datavector_in->get_vector(), 1.0/max_value) );
				
				coutAll << "--- scaling data with max value of x: " << max_value << std::endl;
				coutAll.synchronize();
			}

			LOG_FUNC_STATIC_END
		}
		
		template<>
		void VarX::saveCSV(std::string name_of_file, int T, int xdim, int *xmem_arr, int *K_arr, GeneralVector<PetscVector> *datavector, GeneralVector<PetscVector> *gammavector, GeneralVector<PetscVector> *thetavector){
			LOG_FUNC_STATIC_BEGIN

			Timer timer_saveCSV; 
			timer_saveCSV.restart();
			timer_saveCSV.start();
	
			int nproc = GlobalManager.get_size();
			int my_rank = GlobalManager.get_rank();

			int xmem_max = max_array(nproc, xmem_arr);
			int xmem = xmem_arr[my_rank];
			int K = K_arr[my_rank];
	
			int n,t,k,t_mem,i;
			
			/* compute local sizes */
			int Tlocal; /* local part */
			TRY( VecGetLocalSize(datavector->get_vector(),&Tlocal) );
			Tlocal = Tlocal/(double)xdim;

			/* to manipulate with file */
			std::ofstream myfile;
						
			/* open file to write */
			myfile.open(name_of_file.c_str());

			/* write header to file */
			for(n=0; n<xdim; n++){
				myfile << "x" << n << "_orig,";
			}
			for(k=0; k<K; k++){
				myfile << "gamma" << k << ",";
			}
			for(n=0; n<xdim; n++){
				myfile << "x" << n << "_model";
				if(n+1 < xdim){
					myfile << ",";
				}
			}
			myfile << "\n";

			/* theta */
			int theta_start; /* where in theta start actual coefficients */
			double Ax;
			double x_model_n;
			int blocksize = 1 + xdim*xmem;
			double *theta_arr;
			TRY( VecGetArray(thetavector->get_vector(),&theta_arr) );

			/* gamma */
			double *gamma_arr;
			TRY( VecGetArray(gammavector->get_vector(),&gamma_arr) );

			/* go through processors and write the sequence into local file */
			int t_scatter = 50; /* > xmem */
			int t_in_scatter;
			int is_begin = 0;
			int is_end = min(is_begin + t_scatter,T);

			VecScatter ctx;
			IS scatter_is;
			IS scatter_is_to;
			Vec x_scatter;
				TRY( VecCreateSeq(PETSC_COMM_SELF, t_scatter*xdim, &x_scatter) );
			const double *x_scatter_arr;
			while(is_end <= T && is_begin < is_end){ /* while there is something to scatter */
				/* scatter part of time serie */
				TRY( ISCreateStride(PETSC_COMM_WORLD, (is_end-is_begin)*xdim, is_begin*xdim, 1, &scatter_is) );
				TRY( ISCreateStride(PETSC_COMM_SELF, (is_end-is_begin)*xdim, 0, 1, &scatter_is_to) );

				TRY( VecScatterCreate(datavector->get_vector(),scatter_is, x_scatter,scatter_is_to,&ctx) );
				TRY( VecScatterBegin(ctx,datavector->get_vector(),x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
				TRY( VecScatterEnd(ctx,datavector->get_vector(),x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
				TRY( VecScatterDestroy(&ctx) );

				TRY( ISDestroy(&scatter_is_to) );
				TRY( ISDestroy(&scatter_is) );

				TRY( PetscBarrier(NULL) );
				
				/* write begin of time-serie */
				TRY( VecGetArrayRead(x_scatter, &x_scatter_arr) );
				if(is_begin==0){
					for(t=0;t<xmem_max;t++){
						/* original time_serie */
						for(n=0;n<xdim;n++){
							myfile << x_scatter_arr[t*xdim+n] << ",";
						}
						/* write gamma vectors */
						for(k=0;k<K;k++){
							myfile << "0,";
						}
						/* new time-serie */
						for(n=0;n<xdim;n++){
							myfile << x_scatter_arr[t*xdim+n];
							if(n+1 < xdim){
								myfile << ",";
							}
						}
						myfile << "\n";
					}
				}
				
				/* common times */
				for(t=is_begin+xmem_max;t<is_end;t++){
					t_in_scatter = t - is_begin;
					/* original x */
					for(n=0;n<xdim;n++){
						myfile << x_scatter_arr[t_in_scatter*xdim+n] << ",";
					}
					/* write gamma vectors */
					for(k=0;k<K;k++){
						myfile << gamma_arr[k*(T-xmem)+t-xmem_max] << ",";
					}
					/* compute new time serie from model */
					for(n=0;n<xdim;n++){

						x_model_n = 0;
						/* mu */
						for(k=0;k<K;k++){
							theta_start = k*blocksize*xdim;
							x_model_n += gamma_arr[k*(T-xmem)+t-xmem_max]*theta_arr[theta_start + n*blocksize]; 
						}
				
						/* A */
						for(t_mem = 1; t_mem <= xmem; t_mem++){
							/* add multiplication with A_{t_mem} */
							Ax = 0;
							for(i = 0; i < xdim; i++){
								for(k=0;k<K;k++){
									theta_start = k*blocksize*xdim;
									Ax += gamma_arr[k*(T-xmem)+t-xmem_max]*theta_arr[theta_start + n*blocksize + 1 + (t_mem-1)*xdim + i]*x_scatter_arr[(t_in_scatter-t_mem)*xdim+i]; 
								}
							}
							x_model_n += Ax;
						}

						myfile << x_model_n; //x_scatter_arr[t_in_scatter*xdim+n];
						if(n+1 < xdim){
							myfile << ",";
						}
					}
					myfile << "\n";
				}
				TRY( VecRestoreArrayRead(x_scatter, &x_scatter_arr) );
				
				/* update upper and lower index of scatter_is */
				is_begin = is_begin + t_scatter - xmem_max;
				is_end = min(is_begin + t_scatter,T);

			}
			TRY( VecDestroy(&x_scatter) );

			TRY( VecRestoreArray(gammavector->get_vector(),&gamma_arr) );
			TRY( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

			
			myfile.close();
			TRY(PetscBarrier(NULL));

			/* writing finished */
			timer_saveCSV.stop();
			coutAll <<  " - problem saved to CSV in: " << timer_saveCSV.get_value_sum() << std::endl;
			coutAll.synchronize();

			LOG_FUNC_STATIC_END
		}
#endif
		
		
	} /* end namespace example */
	
} /* end of namespace pascinference */
