/** @file femhat.h
 *  @brief class for reduction and prolongation on fem meshes using hat functions
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEMHAT_H
#define	PASC_FEMHAT_H

#include "general/common/fem.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

/** \class FemHat
 *  \brief Reduction/prolongation between FEM meshes using hat functions.
 *
*/
template<class VectorBase>
class FemHat : public Fem<VectorBase> {
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
		FemHat(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

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
		
		void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		void compute_decomposition_reduced();

};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
#endif



/* ----------------- Fem implementation ------------- */
template <class VectorBase>
FemHat<VectorBase>::FemHat(double fem_reduce) : Fem<VectorBase>(fem_reduce){
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


template<class VectorBase>
FemHat<VectorBase>::FemHat(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce) : Fem<VectorBase>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	this->diff = (decomposition1->get_T() - 1)/(double)(decomposition2->get_T() - 1);

	compute_overlaps();

	LOG_FUNC_END
}

template <class VectorBase>
FemHat<VectorBase>::~FemHat(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template <class VectorBase>
std::string FemHat<VectorBase>::get_name() const {
	return "FEM-HAT";
}

template <class VectorBase>
void FemHat<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(this->is_reduced()) << std::endl;
	output_global <<  " - diff             : " << this->diff << std::endl;
	output_global <<  " - fem_reduce       : " << this->fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	output_global <<  " - overlap" << std::endl;
	output_local <<   "   - left           : " << this->left_overlap << std::endl;
	output_local <<   "     - left_t1_idx  : " << this->left_t1_idx << std::endl;
	output_local <<   "     - left_t2_idx  : " << this->left_t2_idx << std::endl;
	output_local <<   "   - right          : " << this->right_overlap << std::endl;
	output_local <<   "     - right_t1_idx : " << this->right_t1_idx << std::endl;
	output_local <<   "     - right_t2_idx : " << this->right_t2_idx << std::endl;
	output_local.synchronize();
 	
	if(this->decomposition1 == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		this->decomposition1->print(output_global);
		output_global.pop();
	}

	if(this->decomposition2 == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		this->decomposition2->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

template <class VectorBase>
void FemHat<VectorBase>::compute_overlaps() {
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
		this->left_t1_idx = floor(this->diff*(this->decomposition2->get_Tbegin()-1));
	} else {
		this->left_t1_idx = floor(this->diff*(this->decomposition2->get_Tbegin()));
	}
	if(this->right_overlap){
		this->right_t1_idx = floor(this->diff*(this->decomposition2->get_Tend()-1+1));
	} else {
		this->right_t1_idx = floor(this->diff*(this->decomposition2->get_Tend()-1));
	}

	/* compute appropriate indexes in coarse grid */
	if(this->left_overlap){
		this->left_t2_idx = floor((this->decomposition1->get_Tbegin())/this->diff)-1;
	} else {
		this->left_t2_idx = floor(this->decomposition1->get_Tbegin()/this->diff);
	}
	if(this->right_overlap){
		this->right_t2_idx = floor((this->decomposition1->get_Tend()-1)/this->diff)+1;
	} else {
		this->right_t2_idx = floor((this->decomposition1->get_Tend()-1)/this->diff);
	}

	LOG_FUNC_END
}

template<class VectorBase>
void FemHat<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void FemHat<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	//TODO	

	LOG_FUNC_END
}

template<class VectorBase>
void FemHat<VectorBase>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(this->is_reduced()){
		int T_reduced = ceil(this->decomposition1->get_T()*this->fem_reduce);
		
		/* compute new decomposition */
		this->decomposition2 = new Decomposition<VectorBase>(T_reduced, 
				*(this->decomposition1->get_graph()), 
				this->decomposition1->get_K(), 
				this->decomposition1->get_xdim(), 
				this->decomposition1->get_DDT_size(), 
				this->decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->decomposition2 = this->decomposition1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (this->decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (this->decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	this->diff = (this->decomposition1->get_T() - 1)/(double)(this->decomposition2->get_T() - 1);

	compute_overlaps();
	
	LOG_FUNC_END
}


}
} /* end of namespace */

#endif
