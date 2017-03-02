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

/** \class Fem
 *  \brief Reduction/prolongation between FEM meshes using constant functions.
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
		
		double fem_reduce;
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		Fem(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem();

		/** @brief print info about fem
		 * 
		 * @param output where to print
		 */	
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		virtual void reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
		virtual void prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

		//TODO: cannot me used in general FEM!
		double get_diff() const;
		double get_fem_reduce() const;
		virtual std::string get_name() const;

		void set_decomposition_original(Decomposition *decomposition1);
		void set_decomposition_reduced(Decomposition *decomposition2);

		virtual void compute_decomposition_reduced();

		Decomposition* get_decomposition_original() const;
		Decomposition* get_decomposition_reduced() const;

		bool is_reduced() const;
};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
__global__ void kernel_fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
#endif



/* ----------------- Fem implementation ------------- */

Fem::Fem(double fem_reduce){
	LOG_FUNC_BEGIN

	decomposition1 = NULL;
	decomposition2 = NULL;

	this->fem_reduce = fem_reduce;
	this->diff = 0; /* I dont have this information without decompositions */
	
	LOG_FUNC_END
}

Fem::Fem(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce){
	LOG_FUNC_BEGIN

	this->fem_reduce = fem_reduce;

	this->set_decomposition_original(decomposition1);
	this->set_decomposition_reduced(decomposition2);

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_fem_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_fem_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition2->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	LOG_FUNC_END
}

Fem::~Fem(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

std::string Fem::get_name() const {
	return "FEM-SUM";
}

void Fem::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(is_reduced()) << std::endl;
	output_global <<  " - diff             : " << diff << std::endl;
	output_global <<  " - fem_reduce       : " << fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
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
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round(decomposition2->get_Tlocal()*diff), round(decomposition2->get_Tbegin()*diff), 1, &gammak1_sublocal_is) );
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
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round(decomposition2->get_Tlocal()*diff), round(decomposition2->get_Tbegin()*diff), 1, &gammak1_sublocal_is) );
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

double Fem::get_fem_reduce() const {
	return fem_reduce;
}

void Fem::set_decomposition_original(Decomposition *decomposition1) {
	LOG_FUNC_BEGIN

	this->decomposition1 = decomposition1;

	LOG_FUNC_END
}

void Fem::set_decomposition_reduced(Decomposition *decomposition2) {
	LOG_FUNC_BEGIN

	this->decomposition2 = decomposition2;

	LOG_FUNC_END
}

Decomposition* Fem::get_decomposition_original() const {
	return decomposition1;
}

Decomposition* Fem::get_decomposition_reduced() const {
	return decomposition2;
}


void Fem::compute_decomposition_reduced() {
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
		this->set_decomposition_reduced(decomposition1);
	}

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

bool Fem::is_reduced() const {
	bool return_value;

	if(fem_reduce < 1.0) {
		return_value = true;
	} else {
		return_value = false;
	}
	
	return return_value;
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
