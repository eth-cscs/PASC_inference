/** @file fem.h
 *  @brief class for reduction and prolongation on fem meshes
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM_H
#define	PASC_FEM_H

#ifndef USE_PETSCVECTOR
 #error 'FEM is for PETSCVECTOR'
#endif

/* this class is for petscvector */
typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace common {

/** \class FEM manipulation
 *  \brief Reduction/prolongation between FEM meshes.
 *
*/
class Fem {
	protected:
		Decomposition *decomposition1; /**< decomposition of the larger problem */
		Decomposition *decomposition2; /**< decomposition of smaller problem */
		
		#ifdef USE_CUDA
			int blockSize_reduce; /**< block size returned by the launch configurator */
			int minGridSize_reduce; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize_reduce; /**< the actual grid size needed, based on input size */

			int blockSize_prolongate; /**< block size returned by the launch configurator */
			int minGridSize_prolongate; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize_prolongate; /**< the actual grid size needed, based on input size */
		#endif
		
		double diff_reduce;
		double diff_prolongate;
		
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem(Decomposition *decomposition1, Decomposition *decomposition2);

		/** @brief destructor
		*/
		~Fem();
		
		void reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
		void prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

		//TODO: cannot me used in general FEM!
		double get_diff_reduce() const;
		double get_diff_prolongate() const;

};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int local_size, double diff_reduce);
__global__ void kernel_fem_prolongate_data(double *data2, double *data1, int T2, int T1, int local_size, double diff_prolongate);
#endif



/* ----------------- Fem implementation ------------- */

Fem::Fem(Decomposition *decomposition1, Decomposition *decomposition2){
	LOG_FUNC_BEGIN

	this->decomposition1 = decomposition1;
	this->decomposition2 = decomposition2;

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_fem_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_fem_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff_reduce = (decomposition1->get_T())/(double)(decomposition2->get_T());
	diff_prolongate = (decomposition2->get_T()-1.0)/(double)((decomposition1->get_T()-1.0));

	LOG_FUNC_END
}

Fem::~Fem(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

void Fem::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;

	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak1_sublocal_is;
	Vec gammak1_sublocal_Vec;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, (decomposition2->get_Tlocal())*diff_reduce, (decomposition2->get_Tbegin())*diff_reduce, 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
				double mysum = 0.0;
				for(int i=round(t2*diff_reduce); i < round((t2+1)*diff_reduce);i++){
					mysum += gammak1_arr[i];
				}
				gammak2_arr[t2] = mysum;
			}

			TRYCXX( VecRestoreArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_fem_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition2->get_Tlocal(), diff_reduce);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );
		
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );
		TRYCXX( ISDestroy(&gammak1_sublocal_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	LOG_FUNC_END
}

void Fem::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak2_sublocal_is;
	Vec gammak2_sublocal_Vec;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round(decomposition1->get_Tlocal()*diff_prolongate), round(decomposition1->get_Tbegin()*diff_prolongate), 1, &gammak2_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_sublocal_is, &gammak2_sublocal_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_sublocal_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int t1 =0; t1 < decomposition1->get_Tlocal(); t1++){
				int idx = round(t1*diff_prolongate);
				gammak1_arr[t1] = gammak2_arr[idx];
			}

			TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_sublocal_Vec,&gammak2_arr) );
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_sublocal_Vec,&gammak2_arr) );

			kernel_fem_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(gammak2_arr, gammak1_arr, decomposition2->get_T(), decomposition1->get_T(), decomposition1->get_Tlocal(), diff_prolongate);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );
		
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_sublocal_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak2_Vec, gammak2_sublocal_is, &gammak2_sublocal_Vec) );
		TRYCXX( ISDestroy(&gammak2_sublocal_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );
	}

	TRYCXX( VecRestoreArray(gamma1_Vec,&gammak1_arr) );
	TRYCXX( VecRestoreArray(gamma2_Vec,&gammak2_arr) );

	LOG_FUNC_END
}

double Fem::get_diff_reduce() const {
	return diff_reduce;
}

double Fem::get_diff_prolongate() const {
	return diff_prolongate;
}



#ifdef USE_CUDA
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int local_size, double diff_reduce) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < local_size){
		double mysum = 0.0;
		for(int i=round(t2*diff_reduce); i < round((t2+1)*diff_reduce);i++){
			mysum += data1[i];
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_fem_prolongate_data(double *data2, double *data1, int T2, int T1, int local_size, double diff_prolongate) {
	int t1 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t1 < local_size){
		int idx = round(t1*diff_prolongate);
		data1[t1] = data2[idx];
	}
}

#endif




}
} /* end of namespace */

#endif
