/** @file femhat.h
 *  @brief class for reduction and prolongation on fem meshes using hat functions
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEMHAT_H
#define	PASC_FEMHAT_H

#ifndef USE_PETSC
 #error 'FEMHAT is for PETSC'
#endif

#include "common/fem.h"

namespace pascinference {
namespace common {

/** \class FemHat
 *  \brief Reduction/prolongation between FEM meshes using hat functions.
 *
*/
class FemHat : public Fem {
	protected:
		bool left_overlap;		/**< is there overlap to the left side of time axis? */
		bool right_overlap;		/**< is there overlap to the right side of time axis? */
		
		int left_t1_idx;			/**< appropriate left index in fine grid (with overlap) */
		int right_t1_idx;			/**< appropriate right index in fine grid (with overlap) */

		int left_t2_idx;			/**< appropriate left index in coarse grid (with overlap) */
		int right_t2_idx;			/**< appropriate right index in coarse grid (with overlap) */

		void compute_overlaps();

	public:
		/** @brief create FEM mapping between two decompositions
		*/
		FemHat(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		FemHat(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~FemHat();

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		std::string get_name() const;
		
		void reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
		void prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

		void compute_decomposition_reduced();

};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
#endif



/* ----------------- Fem implementation ------------- */
FemHat::FemHat(double fem_reduce) : Fem(fem_reduce){
	LOG_FUNC_BEGIN

	/* some implicit values */
	this->left_overlap = false;
	this->right_overlap = false;

	this->left_t1_idx = -1;
	this->right_t1_idx = -1;
	this->left_t2_idx = -1;
	this->right_t2_idx = -1;

	
	LOG_FUNC_END
}


FemHat::FemHat(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce) : Fem(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff = (decomposition1->get_T() - 1)/(double)(decomposition2->get_T() - 1);

	compute_overlaps();

	LOG_FUNC_END
}

FemHat::~FemHat(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

std::string FemHat::get_name() const {
	return "FEM-HAT";
}

void FemHat::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(is_reduced()) << std::endl;
	output_global <<  " - diff             : " << diff << std::endl;
	output_global <<  " - fem_reduce       : " << fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	output_global <<  " - overlap" << std::endl;
	output_local <<   "   - left           : " << this->left_overlap << std::endl;
	output_local <<   "     - left_t1_idx  : " << this->left_t1_idx << std::endl;
	output_local <<   "     - left_t2_idx  : " << this->left_t2_idx << std::endl;
	output_local <<   "   - right          : " << this->right_overlap << std::endl;
	output_local <<   "     - right_t1_idx : " << this->right_t1_idx << std::endl;
	output_local <<   "     - right_t2_idx : " << this->right_t2_idx << std::endl;
	output_local.synchronize();
 	
	if(decomposition1 == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		decomposition1->print(output_global);
		output_global.pop();
	}

	if(decomposition2 == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		decomposition2->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

void FemHat::compute_overlaps() {
	LOG_FUNC_BEGIN
	
	/* indicator of begin and end overlap */
	if(GlobalManager.get_rank() == 0){
		this->left_overlap = false;
	} else {
		this->left_overlap = true;
	}

	if(GlobalManager.get_rank() == GlobalManager.get_size()-1){
		this->right_overlap = false;
	} else {
		this->right_overlap = true;
	}

	/* compute appropriate indexes in fine grid */
	if(this->left_overlap){
		this->left_t1_idx = floor(this->diff*(decomposition2->get_Tbegin()-1));
	} else {
		this->left_t1_idx = floor(this->diff*(decomposition2->get_Tbegin()));
	}
	if(this->right_overlap){
		this->right_t1_idx = floor(this->diff*(decomposition2->get_Tend()-1+1));
	} else {
		this->right_t1_idx = floor(this->diff*(decomposition2->get_Tend()-1));
	}

	/* compute appropriate indexes in coarse grid */
	if(this->left_overlap){
		this->left_t2_idx = floor((decomposition1->get_Tbegin())/this->diff)-1;
	} else {
		this->left_t2_idx = floor(decomposition1->get_Tbegin()/this->diff);
	}
	if(this->right_overlap){
		this->right_t2_idx = floor((decomposition1->get_Tend()-1)/this->diff)+1;
	} else {
		this->right_t2_idx = floor((decomposition1->get_Tend()-1)/this->diff);
	}

	LOG_FUNC_END
}

void FemHat::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	#ifdef USE_GPU
		TRYCXX( VecCUDACopyToGPU(gamma1_Vec) );
		TRYCXX( VecCUDACopyToGPU(gamma2_Vec) );
	#endif

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
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, right_t1_idx - left_t1_idx, left_t1_idx, 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			int Tbegin2 = decomposition2->get_Tbegin();

			//TODO: OpenMP?
			for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
				double center_t1 = (Tbegin2+t2)*diff;
				double left_t1 = (Tbegin2+t2-1)*diff;
				double right_t1 = (Tbegin2+t2+1)*diff;
				
				int id_counter = floor(left_t1) - left_t1_idx; /* first index in provided local t1 array */

				double phi_value; /* value of basis function */

				/* left part of hat function */
				double mysum = 0.0;
				int t1 = floor(left_t1);

				/* compute linear combination with coefficients given by basis functions */
				while(t1 <= center_t1){
					phi_value = (t1 - left_t1)/(center_t1 - left_t1);
					if(id_counter >= 0){
						mysum += phi_value*gammak1_arr[id_counter];
					}
					t1 += 1;
					id_counter += 1;
				}

				/* right part of hat function */
				while(t1 < right_t1){
					phi_value = (t1 - right_t1)/(center_t1 - right_t1);
					if(id_counter < right_t1_idx - left_t1_idx){
						mysum += phi_value*gammak1_arr[id_counter];
					}
					t1 += 1;
					id_counter += 1;
				}

				gammak2_arr[t2] = mysum;
			}
		
			TRYCXX( VecRestoreArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_femhat_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition1->get_Tbegin(), decomposition2->get_Tbegin(), decomposition1->get_Tlocal(), decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
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

void FemHat::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	#ifdef USE_GPU
		TRYCXX( VecCUDACopyToGPU(gamma1_Vec) );
		TRYCXX( VecCUDACopyToGPU(gamma2_Vec) );
	#endif

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

		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, right_t2_idx - left_t2_idx + 1, left_t2_idx, 1, &gammak2_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_sublocal_is, &gammak2_sublocal_Vec) );

		#ifndef USE_CUDA
			TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_sublocal_Vec,&gammak2_arr) );

			int Tbegin1 = decomposition1->get_Tbegin();

			//TODO: OpenMP?
			for(int t1=0; t1 < decomposition1->get_Tlocal(); t1++){
				int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
				int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

				double t1_left = t2_left_id_orig*diff;
				double t1_right = t2_right_id_orig*diff;

				int t2_left_id = t2_left_id_orig - left_t2_idx;
				int t2_right_id = t2_right_id_orig - left_t2_idx;

				/* value of basis functions */
				double t1_value = 0.0;
				double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
				t1_value += phi_value_left*gammak2_arr[t2_right_id];
				
				double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
				t1_value += phi_value_right*gammak2_arr[t2_left_id];

				gammak1_arr[t1] = t1_value;
			}

			TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_sublocal_Vec,&gammak2_arr) );
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_sublocal_Vec,&gammak2_arr) );

			kernel_femhat_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition1->get_Tbegin(), decomposition2->get_Tbegin(), decomposition1->get_Tlocal(), decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
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

	LOG_FUNC_END
}

void FemHat::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(is_reduced()){
		int T_reduced = ceil(decomposition1->get_T()*fem_reduce);
		
		/* compute new decomposition */
		decomposition2 = new Decomposition(T_reduced, 
				*(decomposition1->get_graph()), 
				decomposition1->get_K(), 
				decomposition1->get_xdim(), 
				decomposition1->get_DDT_size(), 
				decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		decomposition2 = decomposition1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff = (decomposition1->get_T() - 1)/(double)(decomposition2->get_T() - 1);

	compute_overlaps();
	
	LOG_FUNC_END
}



#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double center_t1 = (Tbegin2+t2)*diff;
		double left_t1 = (Tbegin2+t2-1)*diff;
		double right_t1 = (Tbegin2+t2+1)*diff;
				
		int id_counter = floor(left_t1) - left_t1_idx; /* first index in provided local t1 array */

		double phi_value; /* value of basis function */

		/* left part of hat function */
		double mysum = 0.0;
		int t1 = floor(left_t1);

		/* compute linear combination with coefficients given by basis functions */
		while(t1 <= center_t1){
			phi_value = (t1 - left_t1)/(center_t1 - left_t1);
			if(id_counter >= 0){
				mysum += phi_value*data1[id_counter];
			}
			t1 += 1;
			id_counter += 1;
		}

		/* right part of hat function */
		while(t1 < right_t1){
			phi_value = (t1 - right_t1)/(center_t1 - right_t1);
			if(id_counter < right_t1_idx - left_t1_idx){
				mysum += phi_value*data1[id_counter];
			}
			t1 += 1;
			id_counter += 1;
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff) {
	int t1 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t1 < T1local){
		int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
		int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

		double t1_left = t2_left_id_orig*diff;
		double t1_right = t2_right_id_orig*diff;

		int t2_left_id = t2_left_id_orig - left_t2_idx;
		int t2_right_id = t2_right_id_orig - left_t2_idx;

		/* value of basis functions */
		double t1_value = 0.0;
		double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
		t1_value += phi_value_left*data2[t2_right_id];
				
		double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
		t1_value += phi_value_right*data2[t2_left_id];

		data1[t1] = t1_value;
	}
}

#endif




}
} /* end of namespace */

#endif
