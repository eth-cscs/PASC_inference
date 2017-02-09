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
		
		double diff;
		
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
		double get_diff() const;

};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
__global__ void kernel_fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
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
		gridSize_prolongate = (decomposition2->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());

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
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round((decomposition2->get_Tend())*diff)-round((decomposition2->get_Tbegin())*diff)+1, round((decomposition2->get_Tbegin())*diff), 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
				double mysum = 0.0;
				for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
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

			kernel_fem_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition2->get_Tlocal(), diff);
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
	IS gammak1_sublocal_is;
	Vec gammak1_sublocal_Vec;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round((decomposition2->get_Tend())*diff)-round((decomposition2->get_Tbegin())*diff)+1, round((decomposition2->get_Tbegin())*diff), 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		#ifndef USE_CUDA
			TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
				for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
					gammak1_arr[i] = gammak2_arr[t2];
				}
			}

			TRYCXX( VecRestoreArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_fem_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition2->get_Tlocal(), diff);
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

	TRYCXX( VecRestoreArray(gamma1_Vec,&gammak1_arr) );
	TRYCXX( VecRestoreArray(gamma2_Vec,&gammak2_arr) );

	LOG_FUNC_END
}

double Fem::get_diff() const {
	return diff;
}



#ifdef USE_CUDA
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double mysum = 0.0;
		for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
			mysum += data1[i];
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
			data1[i] = data2[t2];
		}
	}
}

#endif




}
} /* end of namespace */

#endif
